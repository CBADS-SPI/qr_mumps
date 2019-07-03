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
!> @file qrm_symbolic.F90
!! This files contains the routine that does the symbolic factorization.
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief This subroutine computes the symbolic QR factorization of a matrix.


subroutine qrm_prune_tree(adata, nth, info)

  use qrm_mem_mod
  use qrm_adata_mod
  use qrm_error_mod
  use qrm_sort_mod
  use qrm_common_mod, protect=>qrm_prune_tree
  implicit none

  type(qrm_adata_type), target    :: adata
  integer                         :: nth
  integer, optional               :: info

  integer                         :: i, p
  integer                         :: c, nlz, leaves, totleaves
  integer                         :: n
  real(kind(1.d0)), allocatable   :: lzero_w(:), proc_w(:)
  real(kind(1.d0))                :: totflops, smallth, rm
  integer, allocatable            :: lzero(:), aux(:)
  integer, pointer                :: parent(:)
  logical                         :: found

  ! error management
  integer                         :: err
  character(len=*), parameter     :: name='qrm_prune_tree'

  err = 0
  
  parent => adata%parent
  
  ! at this point we start going down the tree until we identify a set
  ! of nodes such that the subtrees rooted at them can be scheduled to
  ! threads with a good load balancing. Small nodes (or subtrees) will
  ! be pruned away during the descent.
  
  if(err.eq.0) call qrm_alloc(lzero, adata%nnodes      , err)
  if(err.eq.0) call qrm_alloc(adata%small, adata%nnodes, err)
  if(err.eq.0) call qrm_alloc(lzero_w, adata%nnodes    , err)
  if(err.eq.0) call qrm_alloc(aux, adata%nnodes+2      , err)
  if(err.eq.0) call qrm_alloc(proc_w, nth              , err)
  __QRM_INFO_CHECK(err, name, 'qrm_alloc', 9999)

  smallth = 0.01
10 continue

  totleaves = 0
  adata%small = 0

  totflops = real(adata%weight(adata%nnodes), kind(1.d0))

  !goto 20
  nlz = 0
  ! initialize the l0 layer with the root nodes
  do i=1, adata%nnodes
     if(parent(i) .eq. 0) then
        if(adata%weight(i) .gt. smallth*totflops) then
           nlz = nlz+1
           lzero(nlz)   = i
           lzero_w(nlz) = real(adata%weight(i), kind(1.d0))
        else
           adata%small(i) = 1 ! node is too small; mark it
        end if
     end if
     if(adata%childptr(i+1) .eq. adata%childptr(i)) totleaves = totleaves+1
  end do

  leaves = 0

  ! start the loop 
  godown: do
     if(nth .eq. 1) exit ! shortcut for serial execution
     if(nlz .gt. nth*max(2.d0,(log(real(nth,kind(1.d0)))/log(2.d0))**2)) exit ! exit if already too many nodes in l0
     
     proc_w = 0.d0
     ! sort the nodes in l0 in order of decreasing weight
     call qrm_mergesort(nlz, lzero_w(1:nlz), aux(1:nlz+2), order=-1)
     call qrm_mergeswap(nlz, aux(1:nlz+2), lzero(1:nlz), lzero_w(1:nlz))

     ! map subtrees to threads round-robin 
     do i=1, nlz
        ! find the least loaded proc
        p = minloc(proc_w,1)
        proc_w(p) = proc_w(p) + lzero_w(i)
     end do

     ! all the subtrees have been mapped. Evaluate load balance
     rm = minval(proc_w)/maxval(proc_w)

     if((rm .gt. 0.9) .and. (nlz .ge. 1*nth)) exit ! if balance is higher than 90%, we're happy

     ! if load is not balanced, replace heaviest node with its kids (if any)
     found = .false.
     findn: do
        if(leaves .eq. totleaves) exit godown

        if(leaves .eq. nlz) then
           if(nlz .ge. nth*max(2.d0,(log(real(nth,kind(1.d0)))/log(2.d0))**2)) then 
              exit godown ! all the nodes in l0 are leaves. nothing to do
           else
              smallth = smallth/2.d0
              if(smallth .lt. 0.0001) then
                 exit godown
              else
                 goto 10
              end if
           end if
        end if
        n = lzero(leaves+1) ! n is the node that must be replaced

        ! append children of n 
        do p=adata%childptr(n), adata%childptr(n+1)-1
           c = adata%child(p)
           if(real(adata%weight(c),kind(1.d0)) .gt. smallth*totflops) then
              ! this child is big enough, add it
              found = .true.
              nlz = nlz+1
              lzero  (nlz) = c
              lzero_w(nlz) = real(adata%weight(c),kind(1.d0))
           else
              adata%small(c) = 1 ! node is too smal; mark it
           end if
        end do
        if(found) exit findn ! if at least one child was added then we redo the mapping
        leaves = leaves+1
     end do findn

     ! swap n with last element
     lzero  (leaves+1) = lzero  (nlz)
     lzero_w(leaves+1) = lzero_w(nlz)
     nlz = nlz-1

  end do godown

  ! mark all the children of nodes in l0
  do i=1, nlz
     n = lzero(i)
     do p=adata%childptr(n), adata%childptr(n+1)-1
        c = adata%child(p)
        adata%small(c) = 1
     end do
  end do


20 continue
  
9999 continue
  ! cleanup and return
  call qrm_dealloc(lzero)
  call qrm_dealloc(lzero_w)
  call qrm_dealloc(proc_w)
  call qrm_dealloc(aux)
  if(err.ne.0)  call qrm_dealloc(adata%small)

  if(present(info)) info = err
  return

end subroutine qrm_prune_tree
