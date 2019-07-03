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
!> @file qrm_apply2d.F90
!! This file contains a routine that applies Q or Q' to a set of vectors
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief This function applies Q or Q^T to a single vector
!! 
!! @param[in] qrm_mat   the main qrm data structure after factorization.
!!
!! @param[in] transp    a string saying whether Q or Q^T will be applied. Only the
!!                      first character is important.
!! 
!! @param[in,out] b     a 1d array containing the vector to which Q will 
!!                      be applied. 
!!
subroutine _qrm_apply1d(qrm_mat, transp, b, info)

  use _qrm_spmat_mod
  use qrm_common_mod
  use qrm_string_mod
  use _qrm_solve_mod, protect => _qrm_apply1d
  use _qrm_utils_mod
  implicit none

  type(_qrm_spmat_type), intent(in), target  :: qrm_mat
  _qrm_data, intent(inout)                   :: b(:)
  character(len=*), intent(in)               :: transp
  integer, optional                          :: info
  
  _qrm_data, pointer                         :: pnt(:,:)

  ! Turn the 1D pointer into a 2D one and call the 2D method below
  call _qrm_remap_pnt(b, pnt, size(b,1))

  call _qrm_apply2d(qrm_mat, transp, pnt, info)

  return
   
end subroutine _qrm_apply1d



!> @brief This function applies Q or Q^T to a set of vectors
!! 
!! @param[in] qrm_mat   the main qrm data structure after factorization.
!!
!! @param[in] transp    a string saying whether Q or Q^T will be applied. Only the
!!                      first character is important.
!! 
!! @param[in,out] b     a 2d array containing the vectors to which Q will 
!!                      be applied. 
!!

subroutine _qrm_apply2d(qrm_mat, transp, b, info)

  use _qrm_spmat_mod
  use qrm_common_mod
  use qrm_string_mod
  use _qrm_sdata_mod
  use _qrm_solve_mod, protect => _qrm_apply2d
  use qrm_dscr_mod
#if defined (have_starpu)
  use _qrm_starpu_codelets_mod
  use starpu_f_mod
#endif
  implicit none

  type(_qrm_spmat_type), intent(in), target  :: qrm_mat
  _qrm_data, intent(inout), target           :: b(:,:)
  character(len=*), intent(in)               :: transp
  integer, optional                          :: info

  type(_qrm_rhs_type), allocatable           :: b_rhs(:)
  integer                                    :: nb, nrhs, i, keeph
  integer                                    :: n, inode, node, be, en, in, ret
  type(qrm_dscr_type)                        :: qrm_dscr

  ! error management
  integer                                    :: err
  character(len=*), parameter                :: name='qrm_apply'
  
  err = 0
  
  ! immediately check if the facto was done. Otherwise push an error and return
  if(.not.allocated(qrm_mat%fdata)) then
     err = 14
     call qrm_error_print(err, name)
     goto 9999
  else
     if(.not. qrm_mat%fdata%ok) then
        err = 14
        call qrm_error_print(err, name)
        goto 9999
     end if
  end if

  call qrm_get(qrm_mat, 'qrm_keeph', keeph)
  if(keeph .ne. qrm_yes_) then
     err = 30
     call qrm_error_print(err, name)
     goto 9999
  end if

  ! blocking to deal with multiple rhs
  call qrm_get(qrm_mat, 'qrm_rhsnb', nb)
  nrhs = size(b,2)
  if(nb.le.0) nb = nrhs

  allocate(b_rhs((nrhs-1)/nb+1))

  ! init the descriptor
  call qrm_dscr_init(qrm_dscr)

  do i=1, (nrhs-1)/nb+1
     call _qrm_rhs_init(b_rhs(i), b(:,(i-1)*nb+1:min(i*nb,nrhs)))
     call _qrm_apply_async(qrm_dscr, qrm_mat, transp, b_rhs(i), err)
  end do

  call qrm_barrier(qrm_dscr)
  call qrm_dscr_destroy(qrm_dscr)
  err = merge(err, qrm_dscr%err_status, err.ne.0)
  __QRM_INFO_CHECK(err, name,'qrm_apply_async',9999)
  

9998 continue
  ! cleanup and return
  do i=1, (nrhs-1)/nb+1
     call _qrm_rhs_destroy(b_rhs(i))
  end do
  deallocate(b_rhs)

9999 continue
  if(present(info)) info = err
  return
  
end subroutine _qrm_apply2d



subroutine _qrm_apply_async(qrm_dscr, qrm_mat, transp, b_rhs, info)

  use _qrm_spmat_mod
  use qrm_common_mod
  use qrm_string_mod
  use _qrm_sdata_mod
  use qrm_dscr_mod
  use _qrm_solve_mod, protect => _qrm_apply_async
  use iso_c_binding
#if defined (have_starpu)
  use _qrm_starpu_codelets_mod
  use starpu_f_mod
#endif
  implicit none

  type(qrm_dscr_type)                        :: qrm_dscr
  type(_qrm_spmat_type), intent(in), target  :: qrm_mat
  character(len=*), intent(in)               :: transp
  type(_qrm_rhs_type)                        :: b_rhs
  integer, optional                          :: info

  integer                                    :: i, keeph, n, inode
  integer                                    :: node, be, en, in, ret, nb
  type(qrm_adata_type), pointer              :: adata
  
  
  ! error management
  integer                                    :: err
  character(len=*), parameter                :: name='qrm_apply_async'
  
  err = 0

  if(qrm_dscr%err_status .ne. 0) return
  
  ! simplify 
  adata => qrm_mat%adata

  nb = size(b_rhs%p,2)
  
  ! whether to go top-down or bottom-up
  if(qrm_str_tolower(transp(1:1)) .eq. _qrm_transp) then
     be=1; en=adata%nnodes; in=1
  else
     be=adata%nnodes; en=1; in=-1
  end if


  if(.not.allocated(b_rhs%front_rhs)) then
     allocate(b_rhs%front_rhs(adata%nnodes))
#if defined (have_starpu)
     b_rhs%front_rhs(:)%hdl = c_null_ptr
#endif
  end if

#if defined (have_starpu)
  call starpu_f_matrix_data_register(b_rhs%work%hdl, -1, c_null_ptr, &
                                   & nb, nb, qrm_mat%icntl(qrm_nb_), &
                                   & int(_qrm_sizeof_data,kind=c_size_t))
  do node=1, adata%nnodes
     if(adata%small(node) .ge. 0) then
        if(.not.c_associated(b_rhs%front_rhs(node)%hdl)) then
           call starpu_f_void_data_register(b_rhs%front_rhs(node)%hdl)
        end if
     end if
  end do
#else
    ! allocate workspace
  call qrm_alloc(b_rhs%work%c, nb, qrm_mat%icntl(qrm_nb_), err)
  __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)
#endif

  ! inode is the node index (i.e., the i-th node in the give
  ! postorder) while node is its actual ID
  tree: do inode=be, en, in
     if(adata%small(adata%torder(inode)) .lt. 0) cycle tree
     call _qrm_apply_node_task(qrm_dscr, transp, qrm_mat, inode, b_rhs, err)
     __QRM_INFO_CHECK(err, name,'qrm_apply_node_task',9999)
  end do tree
  
#if defined(have_starpu)
  call starpu_f_data_unregister_submit(b_rhs%work%hdl)
#else
  call qrm_dealloc(b_rhs%work%c)
#endif

  
9999 continue

  if(present(info)) info = err
  return
  
end subroutine _qrm_apply_async


