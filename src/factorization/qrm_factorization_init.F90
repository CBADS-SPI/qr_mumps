!! ##############################################################################################
!!
!! Copyright 2012-2016 CNRS, INPT
!! Copyright 2013-2015 UPS
!!  
!! This file is part of qr_mumps.
!!  
!! qr_mumps is free software: you can redistribute it and/or modify
!! it under the terms of the GNU Lesser General Public License as 
!! published by the Free Software Foundation, either version 3 of 
!! the License, or (at your option) any later version.
!!  
!! qr_mumps is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU Lesser General Public License for more details.
!!  
!! You can find a copy of the GNU Lesser General Public License
!! in the qr_mumps/doc directory.
!!
!! ##############################################################################################


!! ##############################################################################################
!> @file qrm_factorization_init.F90
!! This file contains a subroutine that initializes the factorization
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


#include "qrm_common.h" 
!> @brief This subroutine initializes the data structures needed for
!! the actual factorization.

!> The main task achieved by this routine is the creation of what we
!! call (in mumps terminology) the arrowheads. Basically it builds a
!! list of @link _qrm_fdata_mod::_qrm_front_type @endlink elements
!! (each one corresponding to one front) and associates to each of
!! them the related coefficients of the original matrix in CSR
!! format. This coefficients will be assembled into the front matrix
!! at the moment of its activation (this is done by the @link
!! _qrm_init_front @endlink routine).
!!
!! param[in] qrm_mat the usual blob associated to the problem
!!
subroutine _qrm_factorization_init(qrm_mat, transp, info)
 
  use qrm_mem_mod
  use _qrm_spmat_mod
  use _qrm_fdata_mod
  use _qrm_utils_mod
  use qrm_common_mod
  use qrm_error_mod
  use _qrm_factorization_mod, protect => _qrm_factorization_init
#if defined (have_starpu)
  use starpu_f_mod
#endif
  implicit none
  
  type(_qrm_spmat_type), target  :: qrm_mat
  character, intent(in)          :: transp
  integer, optional              :: info
  
  integer                        :: f, nrowsf, nvalsf, i, j, roff, r, c
  integer                        :: lrow, itmp, k, m, n, cnt, nr, node, nworker
  integer, allocatable           :: icperm(:), irperm(:)
  type(_qrm_front_type), pointer :: front
  type(qrm_adata_type), pointer  :: adata
  type(_qrm_fdata_type), pointer :: fdata
  type(_qrm_spmat_type)          :: mat_tmp
  integer :: ts, te, tr
  integer(kind=qrm_mem_kind)     :: avail_mem, size_arrowheads
  real(kind(1.d0))               :: mem_relax
  
  ! error management
  integer                        :: err
  character(len=*), parameter    :: name='qrm_factorization_init'

  err = 0

  ! Building the arrowhead currently goes through a sorting of the
  ! matric according to cperm and rperm (the first eases the assembly of
  ! the fronts and the second the building of the arrowheads. Maybe it
  ! is possible to spare this expensive operations. 

  if(allocated(qrm_mat%fdata)) then
     call _qrm_fdata_cleanup(qrm_mat%fdata, err)
     __QRM_INFO_CHECK(err, name,'qrm_fdata_cleanup',9999)
  else
     call qrm_fdata_init(qrm_mat%fdata, err)
     __QRM_INFO_CHECK(err, name,'qrm_fdata_init',9999)
  end if  
  allocate(qrm_mat%fdata%front_list(qrm_mat%adata%nnodes))

  adata => qrm_mat%adata
  fdata => qrm_mat%fdata

  if(transp .eq. _qrm_transp) then
     n = qrm_mat%m
     m = qrm_mat%n
  else
     m = qrm_mat%m
     n = qrm_mat%n
  end if
  
  if(err.eq.0) call qrm_alloc(icperm, n, err)
  if(err.eq.0) call qrm_alloc(irperm, m, err)
  __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)

  do i=1, n
     icperm(adata%cperm(i)) = i
  end do
  
  do i=1, m
     irperm(adata%rperm(i)) = i
  end do  

  call _qrm_spmat_alloc(mat_tmp, qrm_mat%nz, m, n, 'coo', info=err)
  __QRM_INFO_CHECK(err, name,'qrm_spmat_alloc',9999)
  
  ! sort the matrix with indices in increasing order wrt rperm and cperm
  if(transp .eq. _qrm_transp) then
     do i=1, qrm_mat%nz
        mat_tmp%jcn(i) = icperm(qrm_mat%irn(i))
        mat_tmp%irn(i) = irperm(qrm_mat%jcn(i))
#if defined(zprec) || defined(cprec)
        mat_tmp%val(i) = conjg(qrm_mat%val(i))
#else
        mat_tmp%val(i) = qrm_mat%val(i)
#endif
     end do
  else
     do i=1, qrm_mat%nz
        mat_tmp%jcn(i) = icperm(qrm_mat%jcn(i))
        mat_tmp%irn(i) = irperm(qrm_mat%irn(i))
        mat_tmp%val(i) = qrm_mat%val(i)
     end do
  end if
  
  call qrm_dealloc(icperm)
  call qrm_dealloc(irperm)

  call _qrm_sort_mat(mat_tmp, rc='row', info=err)
  __QRM_INFO_CHECK(err, name,'qrm_sort_mat',9999)

  size_arrowheads = 0
  fdata%nfronts = adata%nnodes
  
  roff = 1
  i    = 1
  do f = 1, adata%nnodes
     ! for all the rows belonging to this front
     cnt          =  0
     front        => fdata%front_list(f)
     nrowsf       =  adata%stair(f) - roff+1
     front%anrows =  nrowsf
     front%num    =  f

     call qrm_alloc(front%aiptr, max(nrowsf+1,2), err)
     __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)
     size_arrowheads = size_arrowheads + max(nrowsf+1,2)*_qrm_sizeof_i
     front%aiptr(1:2) = 1
     do r = roff, adata%stair(f)
        do while (mat_tmp%irn(i+cnt).eq.r)
           cnt = cnt+1
           if(i+cnt.gt.qrm_mat%nz) exit
        end do
        front%aiptr(r-roff+2) = cnt+1
     end do
     if(err.eq.0) call qrm_alloc(front%ajcn, cnt, err)
     if(err.eq.0) call qrm_alloc(front%aval, cnt, err)
     __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)
     size_arrowheads = size_arrowheads + cnt*_qrm_sizeof_i + cnt*_qrm_sizeof_data
     front%ajcn(1:cnt) = adata%cperm(mat_tmp%jcn(i:i+cnt-1))
     front%aval(1:cnt) = mat_tmp%val(i:i+cnt-1)
     roff = qrm_mat%adata%stair(f)+1
     i = i+cnt

  end do

  call _qrm_spmat_destroy(mat_tmp, all=.true.)
  ! call system_clock(te, tr)
  ! write(*,*)'--->', real(te-ts)/real(tr)

  
  call qrm_get(qrm_mat,'qrm_mem_relax', mem_relax)
  if(mem_relax.lt.0) then
     avail_mem = huge(avail_mem)
  else
     avail_mem = ceiling(real(qrm_mat%gstats(qrm_e_facto_mempeak_),kind(1.d0))*mem_relax,kind=qrm_mem_kind)
  end if
  
  call qrm_facto_mem_init(fdata%ma, avail_mem)
  ! reserve arrowheads memory
  call qrm_facto_mem_get(fdata%ma, size_arrowheads)
  
  
#if defined(have_starpu)
  do node=1, adata%nnodes
     front => fdata%front_list(node)
     if(adata%small(node) .ge. 0) then
        call starpu_f_void_data_register(front%sym_handle)
     else
        front%sym_handle = c_null_ptr
     end if
  end do
  nworker = starpu_f_worker_get_count()
  
  call starpu_f_matrix_data_register(fdata%work%hdl, -1, c_null_ptr, &
       & qrm_mat%icntl(qrm_nb_), &
       & qrm_mat%icntl(qrm_nb_), &
       & qrm_mat%icntl(qrm_nb_), &
       & int(_qrm_sizeof_data,kind=c_size_t))
  
#else
  nworker = 1
  allocate(fdata%work%c(qrm_mat%icntl(qrm_nb_),qrm_mat%icntl(qrm_nb_)))
#endif

  qrm_mat%gstats(qrm_nnz_r_) = 0
  qrm_mat%gstats(qrm_nnz_h_) = 0
  
9999 continue 
  if(err.ne.0) call _qrm_fdata_destroy(qrm_mat%fdata)

  if(present(info)) info = err
  return
 
end subroutine _qrm_factorization_init
