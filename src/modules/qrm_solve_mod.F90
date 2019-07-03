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
!> @file qrm_solve_mod.F90
!! This file contains a module with generic interfaces for all the solve typed routines
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


!> @brief This module contains all the interfaces for the typed routines
!! in the solve phase
module _qrm_solve_mod
  
  interface qrm_apply_async
     subroutine _qrm_apply_async(qrm_dscr, qrm_mat, transp, b_rhs, info)
       use _qrm_spmat_mod
       use _qrm_sdata_mod
       use qrm_dscr_mod
       type(qrm_dscr_type)                  :: qrm_dscr
       type(_qrm_spmat_type), intent(in)    :: qrm_mat
       type(_qrm_rhs_type)                  :: b_rhs
       character(len=*), intent(in)         :: transp
       integer, optional                    :: info
     end subroutine _qrm_apply_async
  end interface qrm_apply_async

  interface qrm_solve_async
     subroutine _qrm_solve_async(qrm_dscr, qrm_mat, transp, b_rhs, x_rhs, info)
       use _qrm_spmat_mod
       use _qrm_sdata_mod
       use qrm_dscr_mod
       type(qrm_dscr_type)                  :: qrm_dscr
       type(_qrm_spmat_type), intent(in)    :: qrm_mat
       type(_qrm_rhs_type)                  :: b_rhs, x_rhs
       character(len=*), intent(in)         :: transp
       integer, optional                    :: info
     end subroutine _qrm_solve_async
  end interface qrm_solve_async

  
  
  !> @brief Generic interface for the
  !! @link                                  ::_qrm_apply @endlink and
  !! @link                                  ::_qrm_apply1d @endlink routines
  !!
  interface qrm_apply
     procedure                              :: _qrm_apply2d, _qrm_apply1d
  end interface qrm_apply


  ! !> @brief Generic interface for the
  ! !! @link                                ::_qrm_apply @endlink and
  ! !! @link                                ::_qrm_apply1d @endlink routines
  ! !!
  interface _qrm_apply
     subroutine _qrm_apply2d(qrm_mat, transp, b, info)
       use _qrm_spmat_mod
       type(_qrm_spmat_type), intent(in)    :: qrm_mat
       _qrm_data, intent(inout)             :: b(:,:)
       character(len=*), intent(in)         :: transp
       integer, optional                    :: info
     end subroutine _qrm_apply2d
     subroutine _qrm_apply1d(qrm_mat, transp, b, info)
       use _qrm_spmat_mod
       type(_qrm_spmat_type), intent(in)    :: qrm_mat
       _qrm_data, intent(inout)             :: b(:)
       character(len=*), intent(in)         :: transp
       integer, optional                    :: info
     end subroutine _qrm_apply1d
  end interface _qrm_apply

  !> @brief Generic interface for the @link ::_qrm_solve_r @endlink routine
  ! interface qrm_solve_r
     ! subroutine _qrm_solve_r(qrm_mat, b, x, rhsb)
       ! use _qrm_spmat_mod
       ! type(_qrm_spmat_type), target      :: qrm_mat
       ! _qrm_data, intent(inout)           :: b(:,:)
       ! _qrm_data, intent(out)             :: x(:,:)
       ! integer, optional                  :: rhsb
     ! end subroutine _qrm_solve_r
  ! end interface qrm_solve_r

  !> @brief Generic interface for the @link ::_qrm_solve_rt @endlink routine
  ! interface qrm_solve_rt
     ! subroutine _qrm_solve_rt(qrm_mat, b, x, rhsb)
       ! use _qrm_spmat_mod
       ! type(_qrm_spmat_type), target      :: qrm_mat
       ! _qrm_data, intent(inout)           :: b(:,:)
       ! _qrm_data, intent(out)             :: x(:,:)
       ! integer, optional                  :: rhsb
     ! end subroutine _qrm_solve_rt
  ! end interface qrm_solve_rt

  !> @brief Generic interface for the
  !! @link                                  ::_qrm_solve @endlink and
  !! @link                                  ::_qrm_solve1d @endlink routines
  !!
  interface qrm_solve
     procedure                              :: _qrm_solve1d, _qrm_solve2d
  end interface qrm_solve

  !> @brief Generic interface for the
  !! @link                                  ::_qrm_solve @endlink and
  !! @link                                  ::_qrm_solve1d @endlink routines
  !!
  interface _qrm_solve
     subroutine _qrm_solve2d(qrm_mat, transp, b, x, info)
       use _qrm_spmat_mod
       type(_qrm_spmat_type), target        :: qrm_mat
       _qrm_data, intent(inout)             :: b(:,:)
       _qrm_data, intent(out)               :: x(:,:)
       character(len=*)                     :: transp
       integer, optional                    :: info
     end subroutine _qrm_solve2d
     subroutine _qrm_solve1d(qrm_mat, transp, b, x, info)
       use _qrm_spmat_mod
       type(_qrm_spmat_type), target        :: qrm_mat
       _qrm_data, intent(inout)             :: b(:)
       _qrm_data, intent(out)               :: x(:)
       character(len=*)                     :: transp
       integer, optional                    :: info
     end subroutine _qrm_solve1d
  end interface _qrm_solve

  !> @brief Generic interface for the @link ::_qrm_solve_sing_front @endlink routine
  ! interface qrm_solve_sing_front
     ! subroutine _qrm_solve_sing_front(qrm_mat, b, x, trans)
       ! use _qrm_spmat_mod
       ! use _qrm_fdata_mod
       ! type(_qrm_spmat_type), target      :: qrm_mat
       ! _qrm_data, intent(inout)           :: b(:,:)
       ! _qrm_data, intent(inout)           :: x(:,:)
       ! character                          :: trans
     ! end subroutine _qrm_solve_sing_front
  ! end interface qrm_solve_sing_front


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Tasks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! interface qrm_apply_asm_task     
     ! subroutine _qrm_apply_asm_task(transp, qrm_mat, front, b, rhsb)
       ! use _qrm_spmat_mod
       ! type(_qrm_spmat_type), intent(in)  :: qrm_mat
       ! _qrm_data, intent(inout)           :: b(:,:)
       ! character(len=*), intent(in)       :: transp
       ! type(_qrm_front_type)              :: front
       ! integer                            :: rhsb
     ! end subroutine _qrm_apply_asm_task
  ! end interface qrm_apply_asm_task

  ! interface qrm_apply_task
     ! subroutine _qrm_apply_task(transp, front, b, work, i)
       ! use _qrm_spmat_mod
       ! _qrm_data, intent(inout)           :: b(:,:)
       ! type(_qrm_bc_type), intent(inout)  :: work
       ! character(len=*), intent(in)       :: transp
       ! type(_qrm_front_type)              :: front
       ! integer                            :: rhsb
     ! end subroutine _qrm_apply_task
  ! end interface qrm_apply_task
  
  ! interface qrm_solve_asm_task     
     ! subroutine _qrm_solve_asm_task(transp, qrm_mat, front, b, x, rhsb)
       ! use _qrm_spmat_mod
       ! type(_qrm_spmat_type), intent(in)  :: qrm_mat
       ! _qrm_data, intent(inout)           :: b(:,:), x(:,:)
       ! character(len=*), intent(in)       :: transp
       ! type(_qrm_front_type)              :: front
       ! integer                            :: rhsb
     ! end subroutine _qrm_solve_asm_task
  ! end interface qrm_solve_asm_task

  ! interface qrm_solve_task
     ! subroutine _qrm_solve_task(transp, front, b, x, i)
       ! use _qrm_spmat_mod
       ! _qrm_data, intent(inout)           :: b(:,:)
       ! _qrm_data, intent(inout)           :: x(:,:)
       ! character(len=*), intent(in)       :: transp
       ! type(_qrm_front_type)              :: front
       ! integer                            :: rhsb
     ! end subroutine _qrm_solve_task
  ! end interface qrm_solve_task

  ! interface qrm_apply_subtree_task
     ! subroutine _qrm_apply_subtree_task(transp, qrm_mat, iroot, b, work, rhsb)
       ! use _qrm_spmat_mod
       ! type(_qrm_spmat_type), intent(in)  :: qrm_mat
       ! _qrm_data, intent(inout)           :: b(:,:)
       ! type(_qrm_bc_type), intent(inout)  :: work
       ! character(len=*), intent(in)       :: transp
       ! integer                            :: iroot, rhsb
     ! end subroutine _qrm_apply_subtree_task
  ! end interface qrm_apply_subtree_task

  ! interface qrm_solve_subtree_task
     ! subroutine _qrm_solve_subtree_task(transp, qrm_mat, iroot, b, x, rhsb)
       ! use _qrm_spmat_mod
       ! type(_qrm_spmat_type), intent(in)  :: qrm_mat
       ! _qrm_data, intent(inout)           :: b(:,:), x(:,:)
       ! character(len=*), intent(in)       :: transp
       ! integer                            :: iroot, rhsb
     ! end subroutine _qrm_solve_subtree_task
  ! end interface qrm_solve_subtree_task

  interface qrm_apply_node_task
     subroutine _qrm_apply_node_task(qrm_dscr, transp, qrm_mat, inode, b, info)
       use _qrm_spmat_mod
       use _qrm_sdata_mod
       use qrm_dscr_mod
       type(qrm_dscr_type)                  :: qrm_dscr
       type(_qrm_spmat_type), target        :: qrm_mat
       type(_qrm_rhs_type)                  :: b
       character(len=*), intent(in)         :: transp
       integer                              :: inode
       integer, optional                    :: info
     end subroutine _qrm_apply_node_task
  end interface qrm_apply_node_task
  
  interface qrm_solve_node_task
     subroutine _qrm_solve_node_task(qrm_dscr, transp, qrm_mat, inode, b, x, info)
       use _qrm_spmat_mod
       use _qrm_sdata_mod
       use qrm_dscr_mod
       type(qrm_dscr_type)                  :: qrm_dscr
       type(_qrm_spmat_type), target        :: qrm_mat
       type(_qrm_rhs_type)                  :: b, x
       character(len=*), intent(in)         :: transp
       integer                              :: inode
       integer, optional                    :: info
     end subroutine _qrm_solve_node_task
  end interface qrm_solve_node_task
  
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Kernels
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  interface qrm_front_qt
     subroutine _qrm_front_qt(front, b, work)
       use _qrm_fdata_mod
       use _qrm_sdata_mod
       type(_qrm_front_type)                :: front
       type(_qrm_rhs_type)                  :: b
       _qrm_data                            :: work(:,:)
     end subroutine _qrm_front_qt
  end interface qrm_front_qt

  interface qrm_assemble_qt
     subroutine _qrm_assemble_qt(qrm_mat, front, b, info)
       use _qrm_spmat_mod
       use _qrm_fdata_mod
       use _qrm_sdata_mod
       type(_qrm_spmat_type)                :: qrm_mat
       type(_qrm_front_type)                :: front
       type(_qrm_rhs_type)                  :: b
       integer, optional                    :: info
     end subroutine _qrm_assemble_qt
  end interface qrm_assemble_qt

  interface qrm_front_r
     subroutine _qrm_front_r(front, b, x)
       use _qrm_fdata_mod
       use _qrm_sdata_mod
       type(_qrm_front_type)                :: front
       type(_qrm_rhs_type)                  :: b, x
     end subroutine _qrm_front_r
  end interface qrm_front_r

  interface qrm_assemble_r
     subroutine _qrm_assemble_r(qrm_mat, front, b, x, info)
       use _qrm_spmat_mod
       use _qrm_sdata_mod
       type(_qrm_spmat_type), target        :: qrm_mat
       type(_qrm_front_type)                :: front
       type(_qrm_rhs_type)                  :: b, x
       integer, optional                    :: info
     end subroutine _qrm_assemble_r
  end interface qrm_assemble_r

  interface qrm_front_q
     subroutine _qrm_front_q(front, b, work)
       use _qrm_fdata_mod
       use _qrm_sdata_mod
       type(_qrm_front_type)                :: front
       type(_qrm_rhs_type)                  :: b
       _qrm_data                            :: work(:,:)
     end subroutine _qrm_front_q
  end interface qrm_front_q

  interface qrm_assemble_q
     subroutine _qrm_assemble_q(qrm_mat, front, b, info)
       use _qrm_spmat_mod
       use _qrm_sdata_mod
       type(_qrm_spmat_type), target        :: qrm_mat
       type(_qrm_front_type)                :: front
       type(_qrm_rhs_type)                  :: b
       integer, optional                    :: info
     end subroutine _qrm_assemble_q
  end interface qrm_assemble_q

  interface qrm_front_rt
     subroutine _qrm_front_rt(front, b, x)
       use _qrm_fdata_mod
       use _qrm_sdata_mod
       type(_qrm_front_type)                :: front
       type(_qrm_rhs_type)                  :: b, x
     end subroutine _qrm_front_rt
  end interface qrm_front_rt

  interface qrm_assemble_rt
     subroutine _qrm_assemble_rt(qrm_mat, front, b, x, info)
       use _qrm_spmat_mod
       use _qrm_sdata_mod
       type(_qrm_spmat_type), target        :: qrm_mat
       type(_qrm_front_type)                :: front
       type(_qrm_rhs_type)                  :: b, x
       integer, optional                    :: info
     end subroutine _qrm_assemble_rt
  end interface qrm_assemble_rt

  interface qrm_apply_subtree
     subroutine _qrm_apply_subtree(transp, qrm_mat, iroot, b, work, info)
       use _qrm_spmat_mod
       use _qrm_sdata_mod
       type(_qrm_spmat_type), intent(in)    :: qrm_mat
       type(_qrm_rhs_type)                  :: b
       _qrm_data, intent(inout)             :: work(:,:)
       character(len=*), intent(in)         :: transp
       integer                              :: iroot
       integer, optional                    :: info
     end subroutine _qrm_apply_subtree
  end interface qrm_apply_subtree

  interface qrm_solve_subtree
     subroutine _qrm_solve_subtree(transp, qrm_mat, iroot, b, x, info)
       use _qrm_spmat_mod
       use _qrm_sdata_mod
       type(_qrm_spmat_type), intent(in)    :: qrm_mat
       type(_qrm_rhs_type)                  :: b, x
       character(len=*), intent(in)         :: transp
       integer                              :: iroot
       integer, optional                    :: info
     end subroutine _qrm_solve_subtree
  end interface qrm_solve_subtree
  
end module _qrm_solve_mod
