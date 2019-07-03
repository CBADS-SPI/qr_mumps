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
!> @file qrm_factorization_core.F90
!! This file holds the core of the QR numerical factorization.
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


#include "qrm_common.h"


!> @brief This is the main factorization routine. It performs the factorization of all
!! the fronts that have been assigned to the node (MPI task). The whole process
!! is OpenMP parallel.

!> This factorization is completely dynamic and results in an out of order execution
!! of computational and symbolic tasks. Each front is cut into block-columns. Computational
!! tasks are scheduled on a per-front or per-block-column basis according to a graph of
!! dependencies; each task is pushed on the queue of the thread who owns the related front.
!! The ownership of a front is fixed at the moment where the front is activated. If a thread
!! runs out of tasks to perform (i.e., its queue is empty) it can steal tasks from other
!! threads (see the @link qrm_task_mod @endlink module for the details).
!!
!! Six different types of tasks can be executed by a thread:
!! - panel    : the panel factorization of a column in a front
!! - update   : the update of a block-column wrt to a panel in a front
!! - activate : activation of a front. this operation consists in computing the full
!!              structure of a front and allocating all the memory areas needed for the
!!              subsequent factorization
!! - assemble : this operation consists in assemblying one block-column in the contribution
!!              block of a front into the father node
!! - clean    : cleans up a front, i.e., it deallocates all the memory areas that are no more
!!              needed
!!
!! @param[in,out] qrm_mat This is the main data structure associated
!!                to a problem. The @link _qrm_factorization_core @endlink
!!                assumes that the analysis has been done and thus that the
!!                qrm_mat%adata field contains meaningful stuff
!!

subroutine _qrm_factorization_core(qrm_dscr, qrm_mat, info)

  use iso_c_binding
  use _qrm_spmat_mod
  use qrm_mem_mod
  use qrm_error_mod
  use qrm_common_mod
  use qrm_dscr_mod
  use _qrm_factorization_mod, protect => _qrm_factorization_core
  
  implicit none

  type(qrm_dscr_type)            :: qrm_dscr
  type(_qrm_spmat_type), target  :: qrm_mat
  integer, optional              :: info

  integer                        :: i, j, p, c, node, fnum, ret
  integer, allocatable           :: work(:)
  type(qrm_adata_type), pointer  :: adata
  type(_qrm_fdata_type), pointer :: fdata
  type(_qrm_front_type), pointer :: front, child
  integer                        :: bc, br, prio

  ! error management
  integer                        :: err
  character(len=*), parameter    :: name='qrm_factorization_core'

  err = 0

  allocate(work(qrm_mat%n))
  
  ! simplify
  adata => qrm_mat%adata
  fdata => qrm_mat%fdata


  front => null()

  ! loop over all nodes
  main: do node=1, adata%nnodes
     ! until a proper error analysis is implemented in starpu we have
     ! to check whether an error was raised by some task in order to
     ! stop further submissions
     if(qrm_dscr%err_status.ne.0) goto 9999

     fnum = adata%torder(node)
     if(adata%small(fnum) .lt. 0) cycle
     prio = (adata%nnodes - node + 1)*4
     front => fdata%front_list(fnum)

     if(adata%small(fnum).gt.0) then
        call qrm_facto_mem_get(fdata%ma, adata%asize(fnum))
        call _qrm_do_subtree_task(qrm_dscr, qrm_mat, fnum, huge(1), err)
        __QRM_INFO_CHECK(err, name,'qrm_do_subtree_task',9999)
     else
        call qrm_facto_mem_get(fdata%ma, adata%asize(front%num))
        call _qrm_activate_front(qrm_mat, front, work, err)
        __QRM_INFO_CHECK(err, name,'qrm_activate_front',9999)
        
        call _qrm_init_task(qrm_dscr, qrm_mat, front%num, huge(1), err)
        __QRM_INFO_CHECK(err, name,'qrm_init_task',9999)
        
        do i = adata%childptr(front%num), adata%childptr(front%num+1)-1
           child => fdata%front_list(adata%child(i))
           if(adata%small(adata%child(i)).eq.0) then
              do br =  (child%npiv)/child%mb + 1, (child%ne-1)/child%mb+1
                 do bc = (child%npiv)/child%nb + 1, child%nc
                    ! if block is below the diagonal it does not belong to the CB; cycle
                    if(min(bc*child%nb,child%n).lt.max((br-1)*child%mb+1,1)) cycle
                    call _qrm_assemble_block_task(qrm_dscr, qrm_mat, child%num, br, bc, prio)
                 end do
              end do
           end if
           ! once the child is assembled, it can be cleaned up
           call _qrm_clean_task(qrm_dscr, qrm_mat, child%num, huge(1), err)
           __QRM_INFO_CHECK(err, name,'qrm_clean_task',9999)
        end do
     
        ! Factorize the front
        call _qrm_factorize_front(qrm_dscr, front, fdata%work, prio)
     end if
     
  end do main

  fnum = adata%torder(adata%nnodes)
  call _qrm_clean_task(qrm_dscr, qrm_mat, fnum, huge(1), err)

9999 continue
  ! cleanup and return
#if defined (have_starpu)
  call starpu_f_data_unregister_submit(fdata%work%hdl)
  fdata%work%hdl = c_null_ptr
#else
  deallocate(fdata%work%c)
#endif
    
  deallocate(work)
  if(present(info)) info = err
  return

end subroutine _qrm_factorization_core


subroutine  _qrm_factorize_front(qrm_dscr, front, work, prio)
  use _qrm_fdata_mod
  use qrm_dscr_mod
  use _qrm_factorization_mod, p=>_qrm_factorize_front
  implicit none
  
  type(qrm_dscr_type)            :: qrm_dscr
  type(_qrm_front_type)          :: front
  integer                        :: prio
  type(_qrm_bc_type)             :: work
  
  integer                        :: lt, ft, lct, i, j, k, l, s, bh

  if(front%ne.eq.0) return
  
  lt = (min(front%m,front%n)-1)/front%nb+1
  
  frontfct: do k=1, lt
     ft = ((k-1)*front%nb)/front%mb+1
     lct = (front%stair(min(k*front%nb, front%n))-1)/front%mb+1
     do i = ft, lct, front%bh
        call _qrm_geqrt_task(qrm_dscr, front, k, i, work, prio+3)
        do j=k+1, front%nc
           call _qrm_gemqrt_task(qrm_dscr, front, k, i, j, work, prio+2)
        end do
        do l=i+1, min(i+front%bh-1,lct)
           call _qrm_tpqrt_task(qrm_dscr, front, k, i, l, 's', work, prio+1)
           do j=k+1, front%nc
              call _qrm_tpmqrt_task(qrm_dscr, front, k, i, l, j, 's', work, prio)
           end do
        end do
     end do
     s = front%bh
     do while (s.le.lct-ft+1)
        do i = ft, lct-s, 2*s
           l = i+s
           if(l.le.lct) then
              call _qrm_tpqrt_task(qrm_dscr, front, k, i, l, 't', work, prio+3)
              do j=k+1, front%nc
                 call _qrm_tpmqrt_task(qrm_dscr, front, k, i, l, j, 't', work, prio)
              end do
           end if
        end do
        s = s*2
     end do
  end do frontfct
  
  return
  
  
end subroutine _qrm_factorize_front
