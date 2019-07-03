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
!> @file qrm_reorder_tree.F90
!! This file contains the routine that computes a reordering of the tree to reduce the search
!! space for task scheduling
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


!> @brief This subroutine reorders the assembly tree in order to reduce
!! the tasks search space.

!> A nice side effect of this is that memory consumption is also reduced. The code here
!! basically follows the idea in:
!!
!! J. W. H. Liu. <em>On the storage requirement in the out-of-core
!!                   multifrontal method for sparse factorization.</em>
!!                   ACM Transactions on Mathematical Software, 12:127–148, 1986.
!!
!! @param[in,out] adata a qrm_adata_type data containing a full
!!                caracterization of the assembly tree and fronts
!!                structure.
!! @TODO Add better description of the algorithm
!!

#include "qrm_common.h"

subroutine qrm_reorder_tree(adata, info)

  use qrm_common_mod, protect=>qrm_reorder_tree
  use qrm_adata_mod
  use qrm_mem_mod
  use qrm_sort_mod
  implicit none

  type(qrm_adata_type)        :: adata
  integer, optional           :: info

  integer, allocatable        :: peaks(:), child_peaks(:), aux(:)
  integer, allocatable        :: child(:), tnch(:)
  integer                     :: root, maxch, i, nl, leaf, f, c, sr, sl, nleaves
  integer                     :: nch, nnodes, node, peak, reord
  integer, parameter          :: qrm_noord_=0, qrm_taskord_=1, qrm_memord_=2

  integer                     :: err
  character(len=*), parameter :: name='qrm_reorder_tree'

  err = 0
  
  reord = qrm_noord_

  if(reord.eq.qrm_noord_) then
     ! no reordertree

  else if(reord.eq.qrm_taskord_) then

     call qrm_alloc(peaks, adata%nnodes, err)
     __QRM_INFO_CHECK(err, name, 'qrm_alloc', 9999)
     
     maxch   = 0
     do node = 1, adata%nnodes
        nch = adata%childptr(node+1)-adata%childptr(node)
        if(nch .gt. maxch) maxch = nch
        peaks(node) = -nch
     end do

     call qrm_alloc(child_peaks, maxch, err)
     __QRM_INFO_CHECK(err, name, 'qrm_alloc', 9999)
     call qrm_alloc(aux, maxch+2, err)
     __QRM_INFO_CHECK(err, name, 'qrm_alloc', 9999)

     nnodes = adata%nnodes
     node   = 1

     main: do leaf=1, adata%nnodes

        nch = adata%childptr(leaf+1)-adata%childptr(leaf)
        if(nch .eq. 0) then
           node = leaf

           ! start going up the subtree
           do
              if(nnodes .eq. 0) exit main
              ! peaks(node).eq.0 only occurs if node is a leaf or
              ! if its peak can be computed
              if(peaks(node) .ne. 0) exit

              peak=0
              ! sort node's children (if any) and compute its peak
              nch = 0
              if(adata%childptr(node) .ne. adata%childptr(node+1)) then
                 do i = adata%childptr(node), adata%childptr(node+1)-1
                    nch = nch+1
                    child_peaks(nch) = peaks(adata%child(i))
                 end do
                 ! sort
                 call qrm_mergesort(nch, child_peaks(1:nch), aux(1:nch+2), order=-1)

                 call qrm_mergeswap(nch, aux(1:nch+2), child_peaks(1:nch), &
                      & adata%child(adata%childptr(node):adata%childptr(node+1)-1))
                 ! compute node's peak
                 do i=1, nch
                    peak = max(peak, i-1+child_peaks(i))
                 end do
              end if

              peaks(node) = max(peak,nch+1)
              nnodes = nnodes -1

              f = adata%parent(node)
              if (f .ne. 0) then
                 peaks(f) = peaks(f)+1
                 if(peaks(f) .eq. 0) then
                    node = f
                 else
                    cycle main
                 end if
              end if
           end do

        end if

     end do main

     call qrm_dealloc(child_peaks, err)
     __QRM_INFO_CHECK(err, name, 'qrm_dealloc', 9999)
     call qrm_dealloc(aux, err)
     __QRM_INFO_CHECK(err, name, 'qrm_dealloc', 9999)
     call qrm_dealloc(peaks, err)
     __QRM_INFO_CHECK(err, name, 'qrm_dealloc', 9999)

  else if (reord.eq.qrm_memord_) then
     call qrm_alloc(child, adata%nnodes, err)
     __QRM_INFO_CHECK(err, name, 'qrm_alloc', 9999)
     child = 0

     node = adata%nnodes
     ! go down the subtree
     
     tree_loop: do
        
        nch=adata%childptr(node+1) - adata%childptr(node)
        do
           if(child(node) .ge. nch) then
              ! if we get here it means that all the children of node
              ! have been already visited. 
              exit
           end if
           
           ! increase the number of visited children of node and keep
           ! going down the tree
           child(node) = child(node)+1
           c = adata%child(adata%childptr(node)+child(node)-1)
           node = c
           cycle tree_loop
        end do
        ! we got the the end of the branch. Now go back up
        node = adata%parent(node)
        if (node .eq. 0) exit
     end do tree_loop
     call qrm_dealloc(child, err)
     __QRM_INFO_CHECK(err, name, 'qrm_dealloc', 9999)
  end if

  
  call qrm_alloc(tnch, adata%nnodes , err)
  __QRM_INFO_CHECK(err, name, 'qrm_alloc', 9999)
  call qrm_alloc(child, adata%nnodes, err)
  __QRM_INFO_CHECK(err, name, 'qrm_alloc', 9999)
  child = 0

  ! compute the list of leaves to start facto. Even if the whole tree
  ! is reordered, only the leaves in l0 are returned; this may have to
  ! be improved. TODO

  tnch    = 0
  nleaves = 0
  ! count the number of non-small children for each node
  do node=1, adata%nnodes
     nch = adata%childptr(node+1)-adata%childptr(node)
     if(nch.eq.0) nleaves=nleaves+1
     do i=adata%childptr(node), adata%childptr(node+1)-1
        c = adata%child(i)
        if(adata%small(c) .eq. 0) tnch(node) = tnch(node)+1
     end do
  end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  sr = 0
  sl = 0
  i = 0
  call qrm_alloc(adata%leaves, nleaves     , err)
  __QRM_INFO_CHECK(err, name, 'qrm_alloc', 9999)
  call qrm_alloc(adata%torder, adata%nnodes, err)
  __QRM_INFO_CHECK(err, name, 'qrm_alloc', 9999)

  ! at this point adata%small(i) = 1 only for nodes that are root of a
  ! sequential subtree. At the end of the loop below adata%small(i)=k
  ! where
  ! - k=0 for all nodes that do not belong to a sequential subtree

  ! - k=l  (l>0) if i is the root of a sequential subtree. In this case
  !              l points to the position whithin adata%torder of the
  !              first leaf of the sequential subtree
  ! - k=-l (l>0) if i is any other node of a sequential subtree. In this case
  !              -l points to the position whithin adata%torder of the
  !              first leaf of the sequential subtree

  nleaves = 0
  node = adata%nnodes
  ! go down the subtree

  nodes_loop: do
     
     if(adata%small(node) .eq. 1) then
        sr = node
     end if
     
     nch=adata%childptr(node+1) - adata%childptr(node)
     do
        if(child(node) .ge. nch) then
           ! if we get here it means that all the children of node
           ! have been already visited. node can be added to the
           ! torder array
           i = i+1
           adata%torder(i) = node
           if(node.eq.sr) then
              ! this means that node is the root of a sequential
              ! subtree. In this case small(node) is set to a
              ! positive value that points to the first leaf of the
              ! sequential subtree
              if(sl.eq.0) sl = i
              adata%small(node) = sl
              ! reset sl and sr to zero
              sl = 0
              sr = 0
           else if (sr.eq.0) then
              ! node does not belong to a sequential subtree
              if((tnch(node) .eq. 0) .and. (adata%small(node).eq.0)) then
                 ! in this case node is a leaf of the pruned tree
                 ! (i.e., the tree without the sequential subtrees)
                    nleaves = nleaves+1
                    adata%leaves(nleaves) = node
                 end if
           else
              if(sl.eq.0) sl = i ! this is the first node of a
              ! sequential subtree that has
              ! no children, i.e., the first
              ! leaf
              adata%small(node) = -sl ! make small(node) point to
              ! the first leaf in the
              ! subtree
           end if
           exit
        end if
        ! increase the number of visited children of node and keep
        ! going down the tree
        child(node) = child(node)+1
        c = adata%child(adata%childptr(node)+child(node)-1)
        node = c
        cycle nodes_loop
     end do
     ! we got the the end of the branch. Now go back up
     node = adata%parent(node)
     if (node .eq. 0) exit

  end do nodes_loop

  adata%nleaves = nleaves
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

9999 continue
  ! cleanup and return
  call qrm_dealloc(child)
  call qrm_dealloc(tnch)
  if(err.ne.0) then
     call qrm_dealloc(peaks       )
     call qrm_dealloc(child_peaks )
     call qrm_dealloc(aux         )
     call qrm_dealloc(adata%leaves)
     call qrm_dealloc(adata%torder)
  end if

  if(present(info)) info = err
  return

end subroutine qrm_reorder_tree




