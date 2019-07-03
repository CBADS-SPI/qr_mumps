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
!> @file qrm_solve2d.F90
!! This file contains a routine that solves for R or R' against multiple vectors
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


! #include "prec.h"
#include "qrm_common.h"



!> @brief This function solves for R or R' against a single vector
!!
!! @param[in] qrm_mat   the main qrm data structure after factorization.
!!
!! @param[in] transp    a string saying whether R or R^T will be solved
!!                      for. Only the first character is important.
!!
!! @param[in]     b     a 1d array containing the RHS vector
!!
!! @param[out]    x     a 1d array containing the solution vector
!!
subroutine _qrm_solve1d(qrm_mat, transp, b, x, info)

  use _qrm_spmat_mod
  use qrm_error_mod
  use _qrm_fdata_mod
  use qrm_string_mod
  use _qrm_utils_mod
  use _qrm_solve_mod, protect => _qrm_solve1d
  use qrm_error_mod
  implicit none

  type(_qrm_spmat_type), target  :: qrm_mat
  _qrm_data, intent(in)          :: b(:)
  _qrm_data, intent(out)         :: x(:)
  character(len=*)               :: transp
  integer, optional              :: info

  _qrm_data, pointer             :: pntb(:,:), pntx(:,:)
  integer                        :: n, i

  n = size(b,1)
  call _qrm_remap_pnt(b, pntb, n)
  n = size(x,1)
  call _qrm_remap_pnt(x, pntx, n)

  call _qrm_solve2d(qrm_mat, transp, pntb, pntx, info)

  return

end subroutine _qrm_solve1d




!> @brief This function solves for R or R' against multiple vectors
!!
!! @param[in] qrm_mat   the main qrm data structure after factorization.
!!
!! @param[in] transp    a string saying whether R or R^T will be solved
!!                      for. Only the first character is important.
!!
!! @param[in]     b     a 2d array containing the RHS vectors
!!
!! @param[out]    x     a 2d array containing the solution vectors
!!
subroutine _qrm_solve2d(qrm_mat, transp, b, x, info)

  use _qrm_spmat_mod
  use qrm_error_mod
  use _qrm_fdata_mod
  use _qrm_sdata_mod
  use qrm_string_mod
  use _qrm_utils_mod
  use qrm_dscr_mod
  use _qrm_solve_mod, protect => _qrm_solve2d
#if defined (have_starpu)
  use _qrm_starpu_codelets_mod ! TODO: add dep in makefile
  use starpu_f_mod
#endif
  implicit none

  type(_qrm_spmat_type), target    :: qrm_mat
  _qrm_data, intent(inout), target :: b(:,:)
  _qrm_data, intent(out), target   :: x(:,:)
  character(len=*)                 :: transp
  integer, optional                :: info

  type(qrm_dscr_type)              :: qrm_dscr
  integer                          :: i, nb, nrhs, be, en, in
  integer                          :: node, inode, ret, keeph
  type(_qrm_front_type), pointer   :: front
  type(_qrm_rhs_type), allocatable :: b_rhs(:), x_rhs(:)


  ! error management
  integer                          :: err
  character(len=*), parameter      :: name='qrm_solve'

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
  if(keeph .lt. 0) then
     err = 30
     call qrm_error_print(err, name)
     goto 9999
  end if

  ! blocking to deal with multiple rhs
  call qrm_get(qrm_mat, 'qrm_rhsnb', nb)
  nrhs = size(b,2)
  if(nb.le.0) nb = nrhs

  allocate(b_rhs((nrhs-1)/nb+1))
  allocate(x_rhs((nrhs-1)/nb+1))

  ! init the descriptor
  call qrm_dscr_init(qrm_dscr)

  do i=1, (nrhs-1)/nb+1
     call _qrm_rhs_init(b_rhs(i), b(:,(i-1)*nb+1:min(i*nb,nrhs)))
     call _qrm_rhs_init(x_rhs(i), x(:,(i-1)*nb+1:min(i*nb,nrhs)))
     call _qrm_solve_async(qrm_dscr, qrm_mat, transp, b_rhs(i), x_rhs(i), err)
  end do

  call qrm_barrier(qrm_dscr)
  call qrm_dscr_destroy(qrm_dscr)
  err = merge(err, qrm_dscr%err_status, err.ne.0)
  __QRM_INFO_CHECK(err, name,'qrm_solve_async',9998)

9998 continue 
  do i=1, (nrhs-1)/nb+1
     call _qrm_rhs_destroy(b_rhs(i))
     call _qrm_rhs_destroy(x_rhs(i))
  end do
  deallocate(b_rhs, x_rhs)

9999 continue 
  if(present(info)) info = err
  return

end subroutine _qrm_solve2d


subroutine _qrm_solve_async(qrm_dscr, qrm_mat, transp, b_rhs, x_rhs, info)

  use _qrm_spmat_mod
  use qrm_error_mod
  use _qrm_fdata_mod
  use _qrm_sdata_mod
  use qrm_string_mod
  use _qrm_utils_mod
  use qrm_dscr_mod
  use _qrm_solve_mod, protect => _qrm_solve_async
#if defined (have_starpu)
  use _qrm_starpu_codelets_mod ! TODO: add dep in makefile
  use starpu_f_mod
#endif
  implicit none

  type(qrm_dscr_type)                        :: qrm_dscr
  type(_qrm_spmat_type), intent(in), target  :: qrm_mat
  character(len=*), intent(in)               :: transp
  type(_qrm_rhs_type)                        :: b_rhs, x_rhs
  integer, optional                          :: info

  integer                                    :: i, nb, nrhs, be, en, in
  integer                                    :: node, inode, ret
  type(qrm_adata_type), pointer              :: adata


  ! error management
  integer                                    :: err
  character(len=*), parameter                :: name='qrm_solve_async'

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

  ! do not reallocate if already done
  if(.not.allocated(x_rhs%front_rhs)) then
     allocate(x_rhs%front_rhs(adata%nnodes))
#if defined (have_starpu)
     x_rhs%front_rhs(:)%hdl = c_null_ptr
#endif
  end if
  if(.not.allocated(b_rhs%front_rhs)) then
     allocate(b_rhs%front_rhs(adata%nnodes))
#if defined (have_starpu)
     b_rhs%front_rhs(:)%hdl = c_null_ptr
#endif
  end if


#if defined (have_starpu)
  do node=1, adata%nnodes
     if(adata%small(node) .ge. 0) then
        if(.not.c_associated(x_rhs%front_rhs(node)%hdl))  call starpu_f_void_data_register(x_rhs%front_rhs(node)%hdl)
        if(.not.c_associated(b_rhs%front_rhs(node)%hdl))  call starpu_f_void_data_register(b_rhs%front_rhs(node)%hdl)
     end if
  end do
#endif

  ! inode is the node index (i.e., the i-th node in the give
  ! postorder) while node is its actual ID
  tree: do inode=be, en, in
     if(adata%small(adata%torder(inode)) .lt. 0) cycle tree
     call _qrm_solve_node_task(qrm_dscr, transp, qrm_mat, inode, b_rhs, x_rhs, err)
     __QRM_INFO_CHECK(err, name,'qrm_solve_node_task',9999)
  end do tree

9999 continue

  if(present(info)) info = err
  return
  
end subroutine _qrm_solve_async
