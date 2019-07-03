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
!> @file qrm_factorization_mod.F90
!! This file contains the @link _qrm_factorization_mod @endlink with  the generic interfaces
!! for all the factorization related routines.
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


!> @brief This module contains all the generic interfaces for the typed routines
!! in the factorization phase
module _qrm_factorization_mod

  !> @brief Generic interface for the @link     ::_qrm_factorization_init
  !> @endlink routine
  interface qrm_factorization_init
     subroutine _qrm_factorization_init(qrm_mat, transp, info)
       use _qrm_spmat_mod
       type(_qrm_spmat_type)                    :: qrm_mat
       character, intent(in)                    :: transp
       integer, optional                        :: info
     end subroutine _qrm_factorization_init
  end interface

  !> @brief Generic interface for the @link     ::_qrm_factorize
  !> @endlink routine
  interface qrm_factorize
     subroutine _qrm_factorize(qrm_mat, transp, info)
       use _qrm_spmat_mod
       type(_qrm_spmat_type)                    :: qrm_mat
       character, optional, intent(in)          :: transp
       integer, optional                        :: info
     end subroutine _qrm_factorize
  end interface

  !> @brief Generic interface for the @link     ::_qrm_factorize_async
  !> @endlink routine
  interface qrm_factorize_async
     subroutine _qrm_factorize_async(qrm_dscr, qrm_mat, transp, info)
       use _qrm_spmat_mod
       use qrm_dscr_mod
       type(_qrm_spmat_type)                    :: qrm_mat
       character, optional, intent(in)          :: transp
       type(qrm_dscr_type)                      :: qrm_dscr
       integer, optional                        :: info
     end subroutine _qrm_factorize_async
  end interface

  !> @brief Generic interface for the @link     ::_qrm_factorization_core
  !> @endlink routine
  interface qrm_factorization_core
     subroutine _qrm_factorization_core(qrm_dscr, qrm_mat, info)
       use _qrm_spmat_mod
       use qrm_dscr_mod
       type(qrm_dscr_type)                      :: qrm_dscr
       type(_qrm_spmat_type)                    :: qrm_mat
       integer, optional                        :: info
     end subroutine _qrm_factorization_core
  end interface

  !> @brief Generic interface for the @link     ::_qrm_init_front
  !> @endlink routine
  interface qrm_init_front
     subroutine _qrm_init_front(qrm_mat, front, info)
       use _qrm_fdata_mod
       use _qrm_spmat_mod
       type(_qrm_spmat_type), target            :: qrm_mat
       type(_qrm_front_type)                    :: front
       integer, optional                        :: info
     end subroutine _qrm_init_front
  end interface

  !> @brief Generic interface for the @link     ::_qrm_factorize_front
  !> @endlink routine
  interface qrm_factorize_front
     subroutine _qrm_factorize_front(qrm_dscr, front, work, prio)
       use _qrm_fdata_mod
       use qrm_dscr_mod
       type(qrm_dscr_type)                      :: qrm_dscr
       type(_qrm_front_type)                    :: front
       integer                                  :: prio
       type(_qrm_bc_type)                       :: work
     end subroutine _qrm_factorize_front
  end interface qrm_factorize_front

  !> @brief Generic interface for the @link     ::_qrm_factorize_front
  !> @endlink routine
  interface qrm_activate_front
     subroutine _qrm_activate_front(qrm_mat, front, work, info)
       use _qrm_fdata_mod
       use _qrm_spmat_mod
       type(_qrm_spmat_type), target            :: qrm_mat
       type(_qrm_front_type)                    :: front
       integer, optional, target                :: work(:)
       integer, optional                        :: info
     end subroutine _qrm_activate_front
  end interface qrm_activate_front


  !> @brief Generic interface for the @link     ::_qrm_do_subtree
  !> @endlink routine
  interface qrm_do_subtree
     subroutine _qrm_do_subtree(qrm_mat, front, flops, info)
       use _qrm_spmat_mod
       use _qrm_fdata_mod
       type(_qrm_spmat_type), target            :: qrm_mat
       type(_qrm_front_type)                    :: front
       real(kind(1.d0))                         :: flops
       integer, optional                        :: info
     end subroutine _qrm_do_subtree
  end interface

  !> @brief Generic interface for the @link     ::_qrm_clean_front
  !> @endlink routine
  interface qrm_clean_front
     subroutine _qrm_clean_front(qrm_mat, front, info)
       use _qrm_spmat_mod
       use _qrm_fdata_mod
       type(_qrm_spmat_type), target            :: qrm_mat
       type(_qrm_front_type)                    :: front
       integer, optional                        :: info
     end subroutine _qrm_clean_front
  end interface


  !> @brief Generic interface for the @link     :: _qrm_assemble
  !> @endlink routine
  interface qrm_assemble_block
     subroutine _qrm_assemble_block(qrm_mat, fnum, br, bc)
       use _qrm_spmat_mod
       type(_qrm_spmat_type), target            :: qrm_mat
       integer                                  :: fnum, br, bc
     end subroutine _qrm_assemble_block
  end interface qrm_assemble_block


  !> @brief Generic interface for the @link     ::_qrm_factorization_finalize_task
  !> @endlink routine
  ! interface qrm_factorization_finalize_task
     ! subroutine _qrm_factorization_finalize_task(qrm_mat)
       ! use _qrm_spmat_mod
       ! type(_qrm_spmat_type)                  :: qrm_mat
     ! end subroutine _qrm_factorization_finalize_task
  ! end interface

  
  interface qrm_do_subtree_task
     subroutine _qrm_do_subtree_task(qrm_dscr, qrm_mat, rootnum, prio, info)
       use iso_c_binding
       use _qrm_spmat_mod
       use qrm_dscr_mod
       type(qrm_dscr_type)                      :: qrm_dscr
       type(_qrm_spmat_type), target            :: qrm_mat
       integer                                  :: rootnum
       integer, optional                        :: prio
       integer, optional                        :: info
     end subroutine _qrm_do_subtree_task
  end interface qrm_do_subtree_task

  ! submit geqrt task in StarPU
  interface qrm_geqrt_task
     subroutine _qrm_geqrt_task(qrm_dscr, front, bk, bi, work, prio)
       use _qrm_spmat_mod
       use qrm_dscr_mod
       type(qrm_dscr_type)                      :: qrm_dscr
       type(_qrm_front_type), target            :: front
       integer                                  :: fnum, bk, bi
       type(_qrm_bc_type)                       :: work
       integer, optional                        :: prio
     end subroutine _qrm_geqrt_task
  end interface qrm_geqrt_task

  ! submit gemqrt task in StarPU
  interface qrm_gemqrt_task
     subroutine _qrm_gemqrt_task(qrm_dscr, front, bk, bi, bj, work, prio)
       use _qrm_spmat_mod
       use qrm_dscr_mod
       type(qrm_dscr_type)                      :: qrm_dscr
       type(_qrm_front_type), target            :: front
       integer                                  :: fnum, bk, bi, bj
       type(_qrm_bc_type)                       :: work
       integer, optional                        :: prio
     end subroutine _qrm_gemqrt_task
  end interface qrm_gemqrt_task

  ! submit tpqrt task in StarPU
  interface qrm_tpqrt_task
     subroutine _qrm_tpqrt_task(qrm_dscr, front, bk, bi, bl, ts, work, prio)
       use _qrm_spmat_mod
       use qrm_dscr_mod
       type(qrm_dscr_type)                      :: qrm_dscr
       type(_qrm_front_type), target            :: front
       integer                                  :: fnum, bk, bi, bl
       character(kind=c_char), value            :: ts
       type(_qrm_bc_type)                       :: work
       integer, optional                        :: prio
     end subroutine _qrm_tpqrt_task
  end interface qrm_tpqrt_task

  ! submit tpmqrt task in StarPU
  interface qrm_tpmqrt_task
     subroutine _qrm_tpmqrt_task(qrm_dscr, front, bk, bi, bl, bj, ts, work, prio)
       use _qrm_spmat_mod
       use qrm_dscr_mod
       type(qrm_dscr_type)                      :: qrm_dscr
       type(_qrm_front_type), target            :: front
       integer                                  :: fnum, bk, bi, bl, bj
       character(kind=c_char), value            :: ts
       type(_qrm_bc_type)                       :: work
       integer, optional                        :: prio
     end subroutine _qrm_tpmqrt_task
  end interface qrm_tpmqrt_task

  ! submit clean task in StarPU
  interface qrm_clean_task
     subroutine _qrm_clean_task(qrm_dscr, qrm_mat, fnum, prio, info)
       use _qrm_spmat_mod
       use qrm_dscr_mod
       type(qrm_dscr_type)                      :: qrm_dscr
       type(_qrm_spmat_type), target            :: qrm_mat
       integer                                  :: fnum
       integer, optional                        :: prio
       integer, optional                        :: info
     end subroutine _qrm_clean_task
  end interface qrm_clean_task

  ! insert the assembly tasks in StarPU
  interface qrm_assemble_block_task
     subroutine _qrm_assemble_block_task(qrm_dscr, qrm_mat, fnum, br, bc, prio)
       use _qrm_spmat_mod
       use qrm_dscr_mod
       type(qrm_dscr_type)                      :: qrm_dscr
       type(_qrm_spmat_type), target            :: qrm_mat
       integer                                  :: br, bc, fnum
       integer, optional                        :: prio
     end subroutine _qrm_assemble_block_task
  end interface qrm_assemble_block_task

  ! insert the init tasks in StarPU
  interface qrm_init_task
     subroutine _qrm_init_task(qrm_dscr, qrm_mat, fnum, prio, info)
       use _qrm_spmat_mod
       use qrm_dscr_mod
       type(qrm_dscr_type)                      :: qrm_dscr
       type(_qrm_spmat_type), target            :: qrm_mat
       integer                                  :: fnum
       integer, optional                        :: prio
       integer, optional                        :: info
     end subroutine _qrm_init_task
  end interface qrm_init_task
  


end module _qrm_factorization_mod
