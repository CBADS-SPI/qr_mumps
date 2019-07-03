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
!> @file qrm_solve_tasks.F90
!! This file contains all the tasks for the solve phase
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################

#include "qrm_common.h"

! subroutine _qrm_apply_asm_task(transp, qrm_mat, front, b, rhsb)
!   use _qrm_spmat_mod
!   use qrm_common_mod
!   use _qrm_solve_mod, protect => _qrm_apply_asm_task
!   implicit none
  
!   type(_qrm_spmat_type), intent(in) :: qrm_mat
!   _qrm_data, intent(inout)          :: b(:,:)
!   character(len=*), intent(in)      :: transp
!   type(_qrm_front_type)             :: front
!   integer                           :: rhsb

! ! #if defined (have_starpu)
!   ! starpu task submission goes here
! ! #else

!   if(transp .eq. _qrm_transp) then
!      call _qrm_assemble_qt(qrm_mat, front, b, rhsb)
!   else
!      call _qrm_assemble_q(qrm_mat, front, b, rhsb)
!   end if
  
! ! #endif

!   return

! end subroutine _qrm_apply_asm_task


! subroutine _qrm_apply_task(transp, front, b, work, rhsb)
!   use _qrm_spmat_mod
!   use qrm_common_mod
!   use _qrm_solve_mod, protect => _qrm_apply_task
!   implicit none
  
!   _qrm_data, intent(inout)          :: b(:,:)
!   type(_qrm_bc_type), intent(inout) :: work
!   character(len=*), intent(in)      :: transp
!   type(_qrm_front_type)             :: front
!   integer                           :: rhsb

!   integer                           :: info
  
! ! #if defined (have_starpu)
!   ! starpu task submission goes here
! ! #else

!   if(transp .eq. _qrm_transp) then
!      call _qrm_front_qt(front, b, work%c, rhsb, info)
!   else
!      call _qrm_front_q(front, b, work%c, rhsb, info)
!   end if

! ! #endif

!   return
  
! end subroutine _qrm_apply_task



! subroutine _qrm_solve_asm_task(transp, qrm_mat, front, b, x, rhsb)
!   use _qrm_spmat_mod
!   use qrm_common_mod
!   use _qrm_solve_mod, protect => _qrm_solve_asm_task
!   implicit none
  
!   type(_qrm_spmat_type), intent(in) :: qrm_mat
!   _qrm_data, intent(inout)          :: b(:,:), x(:,:)
!   character(len=*), intent(in)      :: transp
!   type(_qrm_front_type)             :: front
!   integer                           :: rhsb

! ! #if defined (have_starpu)
!   ! starpu task submission goes here
! ! #else

!   if(transp .eq. _qrm_transp) then
!      call _qrm_assemble_rt(qrm_mat, front, b, x, rhsb)
!   else
!      call _qrm_assemble_r(qrm_mat, front, b, x, rhsb)
!   end if
  
! ! #endif

!   return

! end subroutine _qrm_solve_asm_task


! subroutine _qrm_solve_task(transp, front, b, x, rhsb)
!   use _qrm_spmat_mod
!   use qrm_common_mod
!   use _qrm_solve_mod, protect => _qrm_solve_task
!   implicit none
  
!   _qrm_data, intent(inout)          :: b(:,:)
!   _qrm_data, intent(inout)          :: x(:,:)
!   character(len=*), intent(in)      :: transp
!   type(_qrm_front_type)             :: front
!   integer                           :: rhsb

!   integer                           :: info
  
! ! #if defined (have_starpu)
!   ! starpu task submission goes here
! ! #else

!   if(transp .eq. _qrm_transp) then
!      call _qrm_front_rt(front, b, x, rhsb, info)
!   else
!      call _qrm_front_r(front, b, x, rhsb, info)
!   end if

! ! #endif

!   return
  
! end subroutine _qrm_solve_task


! subroutine _qrm_apply_subtree_task(transp, qrm_mat, iroot, b, work, rhsb)
!   use _qrm_spmat_mod
!   use qrm_common_mod
!   use _qrm_solve_mod, protect => _qrm_apply_subtree_task
!   implicit none
  
!   type(_qrm_spmat_type), intent(in) :: qrm_mat
!   _qrm_data, intent(inout)          :: b(:,:)
!   type(_qrm_bc_type), intent(inout) :: work
!   character(len=*), intent(in)      :: transp
!   integer                           :: iroot, rhsb

! ! #if defined (have_starpu)
!   ! starpu task submission goes here
! ! #else

!      call _qrm_apply_subtree(transp, qrm_mat, iroot, b, work%c, rhsb)
  
! ! #endif

!   return

! end subroutine _qrm_apply_subtree_task


! subroutine _qrm_solve_subtree_task(transp, qrm_mat, iroot, b, x, rhsb)
!   use _qrm_spmat_mod
!   use qrm_common_mod
!   use _qrm_solve_mod, protect => _qrm_solve_subtree_task
!   implicit none
  
!   type(_qrm_spmat_type), intent(in) :: qrm_mat
!   _qrm_data, intent(inout)          :: b(:,:), x(:,:)
!   character(len=*), intent(in)      :: transp
!   integer                           :: iroot, rhsb

! ! #if defined (have_starpu)
!   ! starpu task submission goes here
! ! #else

!      call _qrm_solve_subtree(transp, qrm_mat, iroot, b, x, rhsb)
  
! ! #endif

!   return

! end subroutine _qrm_solve_subtree_task



subroutine _qrm_apply_node_task(qrm_dscr, transp, qrm_mat, inode, b, info)
  ! the Intel compler bugs here when used with vectorization options
  !DIR$ NOOPTIMIZE  
  use _qrm_spmat_mod
  use qrm_common_mod
  use qrm_string_mod
  use _qrm_sdata_mod
  use qrm_dscr_mod
  use _qrm_solve_mod, protect => _qrm_apply_node_task
#if defined (have_starpu)
  use _qrm_starpu_codelets_mod
#endif
  implicit none
  
  type(qrm_dscr_type)               :: qrm_dscr
  type(_qrm_spmat_type), target     :: qrm_mat
  type(_qrm_rhs_type), target       :: b
  integer                           :: inode
  character(len=*), intent(in)      :: transp
  integer, optional                 :: info

  logical                           :: issmall 
  type(_qrm_front_type), pointer    :: front
  integer                           :: node
#if defined (have_starpu)
  type(c_ptr), allocatable          :: chandles(:)
  integer                           :: nc, i, c
  character(kind=c_char)            :: transp_c
#endif

  integer                           :: err
  character(len=*), parameter       :: name='qrm_apply_node_task'

  err = 0
  
  node = qrm_mat%adata%torder(inode)
  
#if defined (have_starpu)

  front => qrm_mat%fdata%front_list(node)

  if(qrm_mat%adata%small(node) .eq. 0) then
     nc = qrm_mat%adata%childptr(node+1)-qrm_mat%adata%childptr(node)
  else
     nc = 0
  end if

  allocate(chandles(nc))
  do i=1, nc 
     c = qrm_mat%adata%child(i-1 + qrm_mat%adata%childptr(node))
     chandles(i) = b%front_rhs(c)%hdl
  end do

  transp_c = transp(1:1)

  call _qrm_insert_apply_task_c(qrm_dscr, c_loc(qrm_mat), &
       & c_loc(b),                                        &
       & front%sym_handle,                                &
       & b%front_rhs(node)%hdl,                           &
       & chandles, b%work%hdl,                            &
       & transp_c, inode, nc, 1)
  
  deallocate(chandles)

#else

  issmall = qrm_mat%adata%small(node) .gt. 0
  
  if(issmall) then
     call _qrm_apply_subtree(transp, qrm_mat, inode, b, b%work%c, err)
     __QRM_INFO_CHECK(err, name, 'qrm_apply_subtree', 9999)
  else
     front => qrm_mat%fdata%front_list(node)
     if(qrm_str_tolower(transp(1:1)) .eq. _qrm_transp) then
        call _qrm_assemble_qt(qrm_mat, front, b, err)
        __QRM_INFO_CHECK(err, name, 'qrm_assemble_qt', 9999)
        call _qrm_front_qt(front, b, b%work%c)
     else
        call _qrm_front_q(front, b, b%work%c)
        call _qrm_assemble_q(qrm_mat, front, b, err)
        __QRM_INFO_CHECK(err, name, 'qrm_assemble_q', 9999)
     end if
  end if

#endif

9999 continue
  if(present(info)) info = err
  return

end subroutine _qrm_apply_node_task


subroutine _qrm_solve_node_task(qrm_dscr, transp, qrm_mat, inode, b, x, info)
  ! the Intel compler bugs here when used with vectorization options
  !DIR$ NOOPTIMIZE  
  use _qrm_spmat_mod
  use qrm_common_mod
  use qrm_string_mod
  use _qrm_sdata_mod
  use qrm_dscr_mod
  use _qrm_solve_mod, protect => _qrm_solve_node_task
#if defined (have_starpu)
  use _qrm_starpu_codelets_mod
#endif
  implicit none
  
  type(qrm_dscr_type)               :: qrm_dscr
  type(_qrm_spmat_type), target     :: qrm_mat
  type(_qrm_rhs_type), target       :: b, x
  character(len=*), intent(in)      :: transp
  integer                           :: inode
  integer, optional                 :: info

  logical                           :: issmall
  type(_qrm_front_type), pointer    :: front
  integer                           :: node

#if defined (have_starpu)
  type(c_ptr), allocatable          :: chandles(:)
  integer                           :: nc, i, c
  character(kind=c_char)            :: transp_c
#endif

  integer                           :: err
  character(len=*), parameter       :: name='qrm_solve_node_task'

  err = 0
  
  node = qrm_mat%adata%torder(inode)

#if defined (have_starpu)

  front => qrm_mat%fdata%front_list(node)

  if(qrm_mat%adata%small(node) .eq. 0) then
     nc = qrm_mat%adata%childptr(node+1)-qrm_mat%adata%childptr(node)
  else
     nc = 0
  end if

  allocate(chandles(nc*2))
  do i=1, nc
     c = qrm_mat%adata%child(i-1 + qrm_mat%adata%childptr(node))
     chandles(i) = x%front_rhs(c)%hdl
  end do
  do i=1, nc
     c = qrm_mat%adata%child(i-1 + qrm_mat%adata%childptr(node))
     chandles(nc+i) = b%front_rhs(c)%hdl
  end do

  transp_c = transp(1:1)
  call _qrm_insert_solve_task_c(qrm_dscr,        &
       & c_loc(qrm_mat), c_loc(b), c_loc(x),     &
       & front%sym_handle,                       &
       & x%front_rhs(node)%hdl,                  &
       & b%front_rhs(node)%hdl,                  &
       & chandles,                               &
       & transp_c, inode, nc , 1)
  
  deallocate(chandles)

#else

  issmall = qrm_mat%adata%small(node) .gt. 0
  
  if(issmall) then
     call _qrm_solve_subtree(transp, qrm_mat, inode, b, x, err)
     __QRM_INFO_CHECK(err, name, 'qrm_solve_subtree', 9999)
  else
     front => qrm_mat%fdata%front_list(node)
     if(qrm_str_tolower(transp(1:1)) .eq. _qrm_transp) then
        call _qrm_assemble_rt(qrm_mat, front, b, x, err)
        __QRM_INFO_CHECK(err, name, 'qrm_assemble_rt', 9999)
        call _qrm_front_rt(front, b, x)
     else
        call _qrm_front_r(front, b, x)
        call _qrm_assemble_r(qrm_mat, front, b, x, err)
        __QRM_INFO_CHECK(err, name, 'qrm_assemble_r', 9999)
     end if
  end if

#endif

9999 continue
  if(present(info)) info = err
  return

end subroutine _qrm_solve_node_task
