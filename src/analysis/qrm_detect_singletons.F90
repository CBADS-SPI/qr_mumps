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
!> @file qrm_detect_singletons.F90
!! this file contains the routine that performs the detection of singletons
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief This subroutine detects singletons in a matrix. 

!> For more deatils please refer to:
!!
!! [1] T. Davis <em>"Multifrontal multithreaded rank-revealing sparse QR factorization"</em>
!!     Accepted for publication on ACM TOMS, 2010
!!
!! The following code largely gathers from the SPQR package singleton detection.
!! Once the singletons are detected, the graph is purged by eliminating all the
!! singleton rows and columns
!!
!! @param[in] graph   the input graph by columns (a qrm_spmat in csc format)
!!
!! @param[out] scol   it is a pointer to integer(:) and on output it
!!                    will return the set of singletons. Need not be
!!                    allocated on entry but on exit it will be of
!!                    size graph%n
!!
!! @param[out] srow   it is a pointer to integer(:) and on output it
!!                    will return the column singletons.  Need not be
!!                    allocated on entry but on exit it will be of
!!                    size graph%n
!!
!! @param[out] mrperm an integer array storing the mapping to the row
!!                    indices of the purged matrix. mrperm(j) = k
!!                    means that column j of the purged graph
!!                    corresponds to column k in the original matrix.
!!                    A row permutation is later computed on the
!!                    purged matrix and, thus, mrperm will be used to
!!                    map this permutation back on the original
!!                    matrix.
!!
!! @param[out] mcperm an integer array storing the mapping to the row
!!                    indices of the purged matrix. mcperm(j) = k
!!                    means that row j of the purged graph corresponds
!!                    to row k in the original matrix.  Because the
!!                    ordering method will be executed on the purged
!!                    graph, this array is necessary to map the
!!                    computed ordering onto the columns of A.
!!
!! @param[out] ncsing the number of column singletons found
!!
!! @param[out] nrsing the number of row singletons found
!!
subroutine _qrm_detect_singletons(graph, scol, srow, mrperm, mcperm, nrsing, ncsing)

  use _qrm_spmat_mod
  use qrm_common_mod
  use qrm_mem_mod

  implicit none

  type(_qrm_spmat_type) :: graph
  integer, allocatable   :: scol(:), srow(:)
  integer, allocatable   :: mcperm(:), mrperm(:)
  integer                :: ncsing, nrsing


  type(_qrm_spmat_type) :: graph_t
  
  integer, allocatable   :: nel(:), isrow(:), tmp(:)
  integer                :: j, row, col, jj, ii, i, curr, ncol_p, nrow_p, nz_p, nz_save, idx
  ! error management
  integer                         :: err_act
  character(len=*), parameter     :: name='qrm_detect_singletons'

  call qrm_err_act_save(err_act)


  ! nel contains the number of nnz in each column
  ! isrow is the row corresponding to a singleton column (singleton row)
  ! squeue is the singletons queue.
  ! srow may not be big enough for internal work so we must allocate isrow
  ! and then copy/compress its content back into srow
  call qrm_realloc(isrow , graph%n)
  call qrm_realloc(mcperm, graph%n)
  call qrm_realloc(mrperm, graph%m)
  call qrm_realloc(scol  , graph%n)
  call qrm_realloc(srow  , graph%m)
  call move_alloc(mcperm, nel)
  __QRM_CHECK_RET(name,'qrm_realloc',9999)


  ! the singleton detection will proceed like this:
  ! a singleton queue will be built with a preliminary
  ! pass over the columns of graph. 
  ! Then, for each singleton in the queue, the corresponding row
  ! is removed (this is equivalent to updating nel on the
  ! columns corresponding to nnzs in the singleton row). Once 
  ! the row is removed, if a new singleton is found, it is added
  ! at the end of the queue.

  ! first pass to fill the queue and compute the degrees
  ncsing = 0
  do j=1, graph%n
     nel(j) = graph%jptr(j+1)-graph%jptr(j)
     if (nel(j) .eq. 0) then
        ! j is a "dead" singleton
        ncsing = ncsing+1
        scol(ncsing) = j
        isrow(ncsing)   = 0
     else if (nel(j) .eq. 1) then
        ! j is a singleton
        ncsing = ncsing+1
        scol(ncsing) = j
        isrow(ncsing)   = graph%irn(graph%jptr(j))
     end if
  end do

  if((ncsing .eq. 0) .or. (ncsing .eq. graph%n)) then
     ! shortcut!
     call qrm_dealloc(isrow)
     __QRM_CHECK_RET(name,'qrm_dealloc',9999)
     call move_alloc(nel, mcperm)
     return
  end if
  
  ! we need both a copy by columns and a copy by rows of the graph
  call _qrm_spmat_convert(graph, graph_t, 'csr', values=.false.)
  __QRM_CHECK_RET(name,'qrm_spmat_convert',9999)
        
  ! loop on all the elements in the queue
  curr = 1
  do
     if(curr .gt. ncsing) exit 
     col = scol(curr)
     row = isrow(curr)
     if ((row .eq. 0)) then
        curr = curr+1
        cycle
     end if
     if(graph_t%iptr(row) .lt.0) then
        ! row was already eliminated so col
        ! has become a dead singleton
        isrow(curr) = 0
        curr = curr+1
        cycle
     end if
     graph_t%iptr(row) = flip(graph_t%iptr(row))
     ! loop over all the elements j in row i
     do jj= unflip(graph_t%iptr(row)), unflip(graph_t%iptr(row+1)) -1
        j = graph_t%jcn(jj)
        ! update the degree
        nel(j) = nel(j)-1
        if (nel(j) .eq. 1) then
           ! new singleton found
           ncsing = ncsing+1
           scol(ncsing) = j
           ! now we must find the corresponding singleton row
           ! we loop over all the elements in the column and look
           ! for that whose row pointer was not flipped in graph_t
           do ii=graph%jptr(j), graph%jptr(j+1)-1
              i = graph%irn(ii)
              if(graph_t%iptr(i) .gt. 0) exit 
           end do
           isrow(ncsing) = i
        end if
     end do
     curr = curr+1
  end do

  ! remove all dead singletons from isrow
  curr=0
  do i=1, ncsing
     if(isrow(i) .ne. 0) then
        curr = curr+1
        srow(curr) = isrow(i)
     end if
  end do

  ! at this point we must prune the graph by removing all the singleton columns
  ! and the corresponding singleton rows.

  ! (non dead) singleton rows are those for which graph_t%iptr is negative
  ! to identify (non dead) singleton columns we must built an inverse permutation array
  ! which (hopefully) will come handy later on
  call move_alloc(nel, mcperm)

  ! purge
  mcperm = 0
  do j = 1,ncsing
     idx = scol(j)
     mcperm(idx) = -j
  end do

  mrperm = 0
  nrow_p = 0
  do i=1, graph_t%m
     if(graph_t%iptr(i) .lt. 0) cycle
     nrow_p = nrow_p+1
     mrperm(i) = nrow_p
  end do

  ncol_p = 0
  nz_p   = 1
  nz_save = 1
  do j=1, graph%n
     ! if mcperm < 0 j is a column singleton. just skip it
     if (mcperm(j) .lt. 0) cycle
     ncol_p = ncol_p+1
     ! this next instruction works because j>=ncol_p
     mcperm(ncol_p) = j
     do ii= graph%jptr(j), graph%jptr(j+1)-1
        i = graph%irn(ii)
        ! go down the column and skip all the elements belonging to singleton rows
        if (graph_t%iptr(i) .lt. 0) cycle
        graph%irn(nz_p) = mrperm(i)
        nz_p = nz_p+1
     end do
     graph%jptr(ncol_p) = nz_save
     nz_save = nz_p
  end do
  graph%jptr(ncol_p+1) = nz_save

#if defined (debug)
  if((ncsing + ncol_p) .ne. graph%n) then
     __QRM_PRNT_DBG('("Inconsistency in singleton detection")')
  end if
#endif

  nrow_p=0
  do i=1, graph_t%m
     if(mrperm(i) .gt. 0) then
        ! this next instruction works because i>=nrow_p
        nrow_p = nrow_p+1
        mrperm(nrow_p) = i
     end if
  end do
  mrperm(nrow_p+1:graph_t%m) = 0


  nrsing   = graph%m-nrow_p
  graph%m  = nrow_p
  graph%n  = ncol_p
  graph%nz = nz_p-1

  call _qrm_spmat_destroy(graph_t, all=.true.)
  call qrm_dealloc(isrow)
  __QRM_CHECK_RET(name,'spmat_destroy/dealloc',9999)

  call qrm_alloc(tmp, ncsing)
  tmp(1:ncsing) = scol(1:ncsing)
  call qrm_dealloc(scol)
  call move_alloc(tmp, scol)

  call qrm_alloc(tmp, nrsing)
  tmp(1:nrsing) = srow(1:nrsing)
  call qrm_dealloc(srow)
  call move_alloc(tmp, srow)


  call qrm_err_act_restore(err_act)
  return

9999 continue ! error management
  call qrm_err_act_restore(err_act)
  if(err_act .eq. qrm_abort_) then
     call qrm_err_check()
  end if

  return
  

contains

  function flip(i)
    integer :: flip
    integer :: i
    flip = -i
    return
  end function flip

  function unflip(i)
    integer :: unflip
    integer :: i
    if (i .lt. 0) then
       unflip = -i
    else
       unflip = i
    end if
    return
  end function unflip


end subroutine _qrm_detect_singletons
