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
!> @file qrm_analyse.F90
!! This file contains the main analysis driver
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief This is the driver routine for the analysis phase.

!> This routine performa a number of symbolic operations in
!! preparation for the numerical factorization:
!! -# it computes the graph of the matrix removing duplicates
!! -# it detects the presence of singletons in the matrix
!! -# it computes a column permutation in order to reduce the fill-in
!! -# it computes the elimination tree
!! -# it postorders the elimination tree
!! -# it computes the rowcount
!! -# it performs an amalgamation of the elimination tree (with fill)
!! -# it computes a row permutation of the matrix in order to have 
!!    a global staircase structure
!! -# merges the singletons into the results of the above operations
!! -# it compresses the results of the above operations to get a size
!!    which is proportional to the number of nodes in the amalgamated tree
!! -# it doea a symbolic factorization which completely characterizes the
!!    structure of fronts etc.
!! -# it computes a tree traversal order which reduces the search space 
!!    for scheduling tasks in the numerical factorization
!!
!! @param[in,out] qrm_mat a qrm_spmat_type data which contains the input matrix.
!!                        On output qrm_mat%adata will contain the results of the 
!!                        analysis phase
!!
!! @param[in] transp      a character saying whether to do analysis on A or A'
!!

subroutine _qrm_analyse(qrm_mat, transp, info)
  use _qrm_spmat_mod
  use qrm_common_mod
  use qrm_error_mod
  use _qrm_utils_mod
  use qrm_string_mod
  use qrm_dscr_mod
  use _qrm_analysis_mod, savesym => _qrm_analyse

  implicit none

  type(_qrm_spmat_type), target   :: qrm_mat
  character, optional             :: transp
  integer, optional               :: info

  type(qrm_dscr_type)             :: qrm_dscr

  integer                         :: err
  character(len=*), parameter     :: name='qrm_analyse'

  err = 0
  
  ! init the descriptor
  call qrm_dscr_init(qrm_dscr)

  ! submit the task(s)
  call _qrm_analyse_async(qrm_dscr, qrm_mat, transp, err)

  ! wait for all the tasks in the dscriptor
  call qrm_barrier(qrm_dscr)
  
  if(present(info)) info = max(err, qrm_dscr%err_status)
  
  ! destroy the descriptor
  call qrm_dscr_destroy(qrm_dscr)
  
  return
end subroutine _qrm_analyse



subroutine _qrm_analyse_async(qrm_dscr, qrm_mat, transp, info)
  use _qrm_spmat_mod
  use qrm_common_mod
  use qrm_error_mod
  use _qrm_utils_mod
  use qrm_string_mod
  use iso_c_binding
  use qrm_dscr_mod
  use _qrm_analysis_mod, savesym => _qrm_analyse_async
#if defined (have_starpu)
  use starpu_f_mod
  use _qrm_starpu_codelets_mod
#endif
  implicit none

  type(qrm_dscr_type)             :: qrm_dscr
  type(_qrm_spmat_type), target   :: qrm_mat
  character, optional             :: transp
  integer, optional               :: info
  
#if defined (have_starpu)
  character(kind=c_char)          :: transp_c
#endif

  ! error management
  integer                         :: err
  character(len=*), parameter     :: name='qrm_analyse_async'

  err = 0
  
  call _qrm_check_spmat(qrm_mat, qrm_analyse_, err)
  __QRM_INFO_CHECK(err, name,'qrm_check_spmat',9999)

  ! make sure all the tasks on this matrix (e.g., the analysis) have completed
  call _qrm_sync_spmat(qrm_mat)
  
  if(allocated(qrm_mat%adata)) then
     call qrm_adata_cleanup(qrm_mat%adata, err)
     __QRM_INFO_CHECK(err, name,'qrm_adata_cleanup',9999)
  else
     call qrm_adata_init(qrm_mat%adata, err)
     __QRM_INFO_CHECK(err, name,'qrm_adata_init',9999)
  end if
  
#if defined(have_starpu)
  if(present(transp)) then
     transp_c = qrm_str_tolower(transp)
  else
     transp_c = 'n'
  end if
  call _qrm_insert_analyse_task_c(qrm_dscr, c_loc(qrm_mat), qrm_mat%adata%hdl, transp_c)
#else
  call _qrm_analysis_core(qrm_mat, qrm_dscr, transp, err)
  __QRM_INFO_CHECK(err, name,'qrm_analysis_core',9999)
#endif

  qrm_mat%adata%ok = .true.
  if(present(info)) info=err
  return

9999 continue
  if(allocated(qrm_mat%adata)) qrm_mat%adata%ok = .false.
  if(present(info)) info=err
  return
end subroutine _qrm_analyse_async



subroutine _qrm_analysis_core(qrm_mat, qrm_dscr, transp, info)

  use _qrm_spmat_mod
  use qrm_common_mod
  use qrm_error_mod
  use _qrm_utils_mod
  use qrm_string_mod
  use qrm_dscr_mod
  use _qrm_analysis_mod, savesym => _qrm_analysis_core

  implicit none

  type(_qrm_spmat_type), target   :: qrm_mat
  type(qrm_dscr_type)             :: qrm_dscr
  character, optional             :: transp
  integer, optional               :: info
  
  type(_qrm_spmat_type), target   :: graph
  integer, allocatable            :: parent(:), rc(:)
  integer, allocatable            :: stair(:), nvar(:)
  integer, allocatable            :: icperm(:), irperm(:)
  integer                         :: ncsing, nrsing
  integer, pointer                :: cperm_p(:), rperm_p(:)
  integer, pointer                :: cperm(:), rperm(:), tmp(:)
  integer                         :: i, nth, bignodes, ilrg, iwrd
  integer                         :: ts, to, te, tp, trc, ta, trp, tc, tsy, tpr, tre, tcm, tr, tg
  character                       :: itransp
  real(kind(1.d0))                :: lrg, wrd, rtmp
  ! error management
  integer                         :: err, err2
  character(len=*), parameter     :: name='qrm_analysis_core'

  err = 0; err2 = 0
  
  __QRM_PRNT_DBG('("Entering the analysis driver")')

  ncsing=0
  nrsing=0

  nullify(cperm_p,rperm_p,cperm,rperm,tmp)

  if(present(transp)) then
     itransp = qrm_str_tolower(transp)
  else
     itransp = 'n'
  end if
  call system_clock(ts)
  call _qrm_compute_graph(qrm_mat, graph, itransp, err)
  __QRM_INFO_CHECK(err, name,'qrm_compute_graph',9999)

  call qrm_adata_init(graph%adata, register=.false.)
  
  ! the cperm array will contain
  call qrm_alloc(graph%adata%cperm, graph%n, err)
  __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)
  call qrm_alloc(graph%adata%rperm, graph%m, err)
  __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)

  cperm => graph%adata%cperm
  rperm => graph%adata%rperm

  call system_clock(tg)
  ! Time to go for the ordering
  call _qrm_do_ordering(graph, cperm, qrm_mat%cperm_in, err)
  __QRM_INFO_CHECK(err, name,'qrm_do_ordering',9999)
  call system_clock(to)
  
  ! build the elimination tree
  call qrm_alloc(parent, graph%n, err)
  __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)

  call _qrm_elim_tree(graph, cperm, parent, err)
  __QRM_INFO_CHECK(err, name,'qrm_elim_tree',9999)
  call system_clock(te)

  ! compute a postorder traversal of the tree
  call qrm_postorder(parent, graph%n, cperm, info=err)
  __QRM_INFO_CHECK(err, name,'qrm_postorder',9999)
  call system_clock(tp)

  call qrm_alloc(rc, graph%n, err)
  __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)

  ! do the symbolic facto
  call _qrm_rowcount(graph, parent, cperm, rc, err)
  __QRM_INFO_CHECK(err, name,'qrm_rowcount',9999)
  call system_clock(trc)

  call qrm_alloc(nvar, graph%n, err)
  __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)

  ! amalgamate the tree
  call qrm_amalg_tree(graph%n, parent, rc, cperm, nvar, &
       &              qrm_mat%icntl(qrm_minamalg_), qrm_mat%rcntl(1), err)
  __QRM_INFO_CHECK(err, name,'qrm_amalg_tree',9999)
  call system_clock(ta)

  call qrm_alloc(stair, graph%n, err)
  __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)

  ! Compute the row permutation to put the matrix in staircase form
  call _qrm_rowperm(graph, cperm, rperm, nvar, stair, err)
  __QRM_INFO_CHECK(err, name,'qrm_rowperm',9999)
  call system_clock(trp)
  
  ! Rework all the data and compress it to nnodes size instead of n
  call qrm_compress_data(graph%adata, cperm, parent, rc, stair, graph%n, err)
  __QRM_INFO_CHECK(err, name,'qrm_compress_data',9999)
  call system_clock(tc)

  ! Do the symbolic facto: estimate fronts size, flops, parallelism etc.
  call _qrm_symbolic(graph, err)
  __QRM_INFO_CHECK(err, name,'qrm_symbolic',9999)
  call system_clock(tsy)
  
  nth = qrm_get_num_threads(qrm_dscr)

  ! find a layer such that all the subtrees rooted there are treated sequentially
  call qrm_prune_tree(graph%adata, nth, err)
  __QRM_INFO_CHECK(err, name,'qrm_prune_tree',9999)
  call system_clock(tpr)

  ! reorder the tree in order to minimize something
  call qrm_reorder_tree(graph%adata, err)
  __QRM_INFO_CHECK(err, name,'qrm_reorder_tree',9999)
  call system_clock(tre)

  call qrm_adata_move(graph%adata, qrm_mat%adata)
  qrm_mat%gstats = graph%gstats

  ! estimate the memory consumption  
  call qrm_compute_memory(qrm_mat, itransp)
  call system_clock(tcm, tr)
  
  ! compute the number of rows in each front and some global stats
  ! call qrm_prnt_array(qrm_mat%adata%cperm(1:qrm_mat%n), 'cp')
  ! call qrm_prnt_array(qrm_mat%adata%cp_ptr(1:qrm_mat%adata%nnodes+1), 'pt')
  ! call qrm_prnt_array(qrm_mat%adata%rperm(1:qrm_mat%m), 'rp')
  qrm_mat%adata%ok = .true.
  bignodes = 0
  lrg      = 0.d0
  wrd      = 0.d0
  
  do i=1, qrm_mat%adata%nnodes
     if(qrm_mat%adata%small(i).eq.0) bignodes = bignodes+1
     if(min(qrm_mat%adata%rc(i),qrm_mat%adata%nfrows(i)).le.0) cycle
     rtmp = real(qrm_mat%adata%rc(i),kind(1.d0))*real(qrm_mat%adata%nfrows(i),kind(1.d0))
     if(rtmp .gt. lrg) then
        lrg  = rtmp
        ilrg = i
     end if
     rtmp = real(qrm_mat%adata%rc(i),kind(1.d0))/real(qrm_mat%adata%nfrows(i),kind(1.d0))
     rtmp = max(rtmp,1.d0/rtmp)
     if(rtmp.gt.wrd) then
        wrd  = rtmp
        iwrd = i
     end if
  end do
  
  __QRM_PRNT_DBG('("  Total ana time                      : ",f7.4)')real(tcm-ts)/real(tr)
  __QRM_PRNT_DBG('("     Graph    time                    : ",f7.4)')real(tg-ts)/real(tr)
  __QRM_PRNT_DBG('("     Ordering time                    : ",f7.4)')real(to-tg)/real(tr)
  __QRM_PRNT_DBG('("     Elimtree time                    : ",f7.4)')real(te-to)/real(tr)
  __QRM_PRNT_DBG('("     Postorde time                    : ",f7.4)')real(tp-te)/real(tr)
  __QRM_PRNT_DBG('("     Rowcount time                    : ",f7.4)')real(trc-tp)/real(tr)
  __QRM_PRNT_DBG('("     Amalgama time                    : ",f7.4)')real(ta-trc)/real(tr)
  __QRM_PRNT_DBG('("     Rowperm  time                    : ",f7.4)')real(trp-ta)/real(tr)
  __QRM_PRNT_DBG('("     Compress time                    : ",f7.4)')real(tc-trp)/real(tr)
  __QRM_PRNT_DBG('("     Symbolic time                    : ",f7.4)')real(tsy-tc)/real(tr)
  __QRM_PRNT_DBG('("     Prune    time                    : ",f7.4)')real(tpr-tsy)/real(tr)
  __QRM_PRNT_DBG('("     Reorder  time                    : ",f7.4)')real(tre-tpr)/real(tr)
  __QRM_PRNT_DBG('("     Compmemo time                    : ",f7.4)')real(tcm-tre)/real(tr)
  __QRM_PRNT_DBG('("  # nodes in the tree                 : ",i10)')qrm_mat%adata%nnodes
  __QRM_PRNT_DBG('("  # big nodes                         : ",i10)')bignodes
  __QRM_PRNT_DBG('("  # nleaves                           : ",i20)')qrm_mat%adata%nleaves  
  __QRM_PRNT_DBG('("  Largest  node (m,n)                 : ",i6,2x,i6)')qrm_mat%adata%nfrows(ilrg), &
       & qrm_mat%adata%rc(ilrg)
  __QRM_PRNT_DBG('("  Weirdest node (m,n)                 : ",i6,2x,i6)')qrm_mat%adata%nfrows(iwrd), &
       & qrm_mat%adata%rc(iwrd)
  __QRM_PRNT_DBG('("  Estimated nonzeroes in R            : ",i20)')qrm_mat%gstats(qrm_e_nnz_r_)
  __QRM_PRNT_DBG('("  Estimated nonzeroes in H            : ",i20)')qrm_mat%gstats(qrm_e_nnz_h_)
  __QRM_PRNT_DBG('("  Estimated total flops at facto      : ",i20)')qrm_mat%gstats(qrm_e_facto_flops_)
  __QRM_PRNT_DBG('("  Estimated memory peak at facto (MB) : ",f20.3)')&
       & real(qrm_mat%gstats(qrm_e_facto_mempeak_))/1e6
  ! call qrm_print_tree('atree.dot',qrm_mat%adata,small=.false.)

9999 continue
  qrm_dscr%err_status = err
  ! cleanup and return
  if(err2.eq.0) call qrm_dealloc(parent, err2)
  if(err2.eq.0) call qrm_dealloc(rc,     err2)
  if(err2.eq.0) call qrm_dealloc(nvar,   err2)
  if(err2.eq.0) call qrm_dealloc(stair,  err2)
  if(err2.eq.0) call qrm_spmat_destroy(graph, all=.true., info=err2)

  if(present(info)) info = merge(err, err2, err.ne.0)
  
  return

end subroutine _qrm_analysis_core


