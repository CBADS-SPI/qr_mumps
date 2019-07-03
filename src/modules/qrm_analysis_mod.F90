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
!> @file qrm_analysis_mod.F90
!! This file contains a module with all the generic interfaces for the typed analysis routines.
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


!> @brief This module contains the generic interfaces for all the analysis routines.
module _qrm_analysis_mod

  !> @brief Generic interface for the @link ::_qrm_analyse @endlink routine
  interface qrm_analyse
     subroutine _qrm_analyse(qrm_mat, transp, info)
       use _qrm_spmat_mod
       type(_qrm_spmat_type)                :: qrm_mat
       character, optional, intent(in)      :: transp
       integer, optional                    :: info
     end subroutine _qrm_analyse
  end interface qrm_analyse

  !> @brief Generic interface for the @link ::_qrm_analyse @endlink routine
  interface qrm_analyse_async
     subroutine _qrm_analyse_async(qrm_dscr, qrm_mat, transp, info)
       use qrm_dscr_mod
       use _qrm_spmat_mod
       type(qrm_dscr_type)                  :: qrm_dscr
       type(_qrm_spmat_type)                :: qrm_mat
       character, optional, intent(in)      :: transp
       integer, optional                    :: info
     end subroutine _qrm_analyse_async
  end interface qrm_analyse_async

  !> @brief Generic interface for the @link ::_qrm_analyse @endlink routine
  interface qrm_analysis_core
     subroutine _qrm_analysis_core(qrm_mat, qrm_dscr, transp, info)
       use qrm_dscr_mod
       use _qrm_spmat_mod
       type(_qrm_spmat_type)                :: qrm_mat
       type(qrm_dscr_type)                  :: qrm_dscr
       character, optional, intent(in)      :: transp
       integer, optional                    :: info
     end subroutine _qrm_analysis_core
  end interface qrm_analysis_core

  !> @brief Generic interface for the @link ::_qrm_compute_graph @endlink routine
  interface qrm_compute_graph
     subroutine _qrm_compute_graph(qrm_mat, graph, transp, info)
       use _qrm_spmat_mod
       type(_qrm_spmat_type)                :: qrm_mat
       type(_qrm_spmat_type), intent(out)   :: graph
       character                            :: transp
       integer, optional                    :: info
     end subroutine _qrm_compute_graph
  end interface

  !> @brief Generic interface for the @link ::_qrm_detect_singletons @endlink routine
  interface qrm_detect_singletons
     subroutine _qrm_detect_singletons(graph, scol, srow, mrperm, mcperm, nrsing, ncsing)
       use _qrm_spmat_mod
       type(_qrm_spmat_type), intent(in)    :: graph
       integer, allocatable                 :: scol(:), srow(:)
       integer, allocatable                 :: mcperm(:), mrperm(:)
       integer                              :: nrsing, ncsing
     end subroutine _qrm_detect_singletons
  end interface

  !> @brief Generic interface for the @link ::_qrm_do_colamd @endlink routine
  interface qrm_do_colamd
     subroutine _qrm_do_colamd(graph, cperm, info)
       use _qrm_spmat_mod
       type(_qrm_spmat_type)                :: graph
       integer                              :: cperm(:)
       integer, optional                    :: info
     end subroutine _qrm_do_colamd
  end interface

  !> @brief Generic interface for the @link ::_qrm_do_metis @endlink routine
  interface qrm_do_metis
     subroutine _qrm_do_metis(graph, cperm, info)
       use _qrm_spmat_mod
       type(_qrm_spmat_type)                :: graph
       integer                              :: cperm(:)
       integer, optional                    :: info
     end subroutine _qrm_do_metis
  end interface

  !> @brief Generic interface for the @link ::_qrm_do_scotch @endlink routine
  interface qrm_do_scotch
     subroutine _qrm_do_scotch(graph, cperm, info)
       use _qrm_spmat_mod
       type(_qrm_spmat_type)                :: graph
       integer                              :: cperm(:)
       integer, optional                    :: info
     end subroutine _qrm_do_scotch
  end interface

  !> @brief Generic interface for the @link ::_qrm_ata_graph @endlink routine
  interface qrm_ata_graph
     subroutine _qrm_ata_graph(g_csc, ata_graph, info)
       use _qrm_spmat_mod
       type(_qrm_spmat_type), intent(in)    :: g_csc
       type(_qrm_spmat_type), intent(out)   :: ata_graph
       integer, optional                    :: info
     end subroutine _qrm_ata_graph
  end interface

  !> @brief Generic interface for the @link ::_qrm_do_ordering @endlink routine
  interface qrm_do_ordering
     subroutine _qrm_do_ordering(graph, cperm, cperm_in, info)
       use _qrm_spmat_mod
       type(_qrm_spmat_type)                :: graph
       integer                              :: cperm(:)
       integer, pointer                     :: cperm_in(:)
       integer, optional                    :: info
     end subroutine _qrm_do_ordering
  end interface

  !> @brief Generic interface for the @link ::_qrm_elim_tree @endlink routine
  interface qrm_elim_tree
     subroutine _qrm_elim_tree(graph, cperm, parent, info)
       use _qrm_spmat_mod
       type(_qrm_spmat_type), intent(in)    :: graph
       integer, intent(in)                  :: cperm(:)
       integer, allocatable                 :: parent(:)
       integer, optional                    :: info
     end subroutine _qrm_elim_tree
  end interface

  !> @brief Generic interface for the @link ::_qrm_rowcount @endlink routine
  interface qrm_rowcount
     subroutine _qrm_rowcount(graph, parent, porder, rc, info)
       use _qrm_spmat_mod
       type(_qrm_spmat_type)                :: graph
       integer                              :: parent(:), porder(:), rc(:)
       integer, optional                    :: info
     end subroutine _qrm_rowcount
  end interface

  !> @brief Generic interface for the @link ::_qrm_symbolic @endlink routine
  interface qrm_symbolic
     subroutine _qrm_symbolic(graph, info)
       use _qrm_spmat_mod
       type(_qrm_spmat_type), intent(inout) :: graph
       integer, optional                    :: info
     end subroutine _qrm_symbolic
  end interface
  
  !> @brief Generic interface for the @link ::_qrm_rowperm @endlink routine
  interface qrm_rowperm
     subroutine _qrm_rowperm(qrm_mat, cperm, rperm, nvar, stair, info)
       use _qrm_spmat_mod
       type(_qrm_spmat_type)                :: qrm_mat
       integer                              :: cperm(:), rperm(:), nvar(:), stair(:)
       integer, optional                    :: info
     end subroutine _qrm_rowperm
  end interface

  !> @brief Generic interface for the @link ::_qrm_attach_singletons @endlink routine
  interface qrm_attach_singletons
     subroutine _qrm_attach_singletons(qrm_mat, scol, srow, mrperm, mcperm, nrsing, ncsing)
       use _qrm_spmat_mod
       type(_qrm_spmat_type):: qrm_mat
       integer                              :: scol(:), srow(:), mcperm(:), mrperm(:)
       integer                              :: nrsing, ncsing
     end subroutine _qrm_attach_singletons
  end interface

  interface qrm_compute_memory
     subroutine _qrm_compute_memory(qrm_mat, transp, info)
       use _qrm_spmat_mod
       type(_qrm_spmat_type), target        :: qrm_mat
       character                            :: transp
       integer, optional                    :: info
     end subroutine _qrm_compute_memory
  end interface qrm_compute_memory

  
end module _qrm_analysis_mod
   
