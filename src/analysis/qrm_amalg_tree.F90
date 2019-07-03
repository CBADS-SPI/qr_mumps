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
!> @file qrm_amalg_tree.F90
!! This file contains the @link qrm_amalg_tree @endlink routine which performs the amalgamation
!! of the elimination tree.
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief This subroutine performs amalgamation on the assembly tree.

!> Amalgamation is done in order to have a better exploitation of
!! level-3 blas routines at the price of some more oprations due to
!! added fill-in. Any node can be amalgamated to its parent depending
!! on the amount of fill-in introduced in R. For the moment the number
!! of flops as well as the fill-on introduced in H are not taken into
!! account because we are not capable of predicting the structure of
!! H.
!!
!! @param[in] n            the size of the tree
!! 
!! @param[in,out] parent   both in input and in out this array will contain the assembly tree
!!                         The meaning of parent on output is:
!!                         - parent(i) = j>0: i is the principal variable of a node and j is the principal
!!                                    variable of its father node
!!                         - parent(i) = j=0: i is a principal variable of a root node
!!                         - parent(i) = j<0: i is a subordinate variable inside a node whose
!!                                    principal variable is j.
!!
!! @param[in,out] rowcount an integer array containing the rowcount for thr R factor. This corresponds
!!                         to the front column-size
!!
!! @param[in,out] porder   the tree postorder. Everything below strongly relies on the fact that in porder
!!                         a principal variable comes always before its subordinates
!!
!! @param[out]    nvar     nvar(i) will be set to the number of variables in the node whose principal
!!                         variable is i
!!                         
!! @param[in] min_var      supernodes containing less than min_var, will be amalgamated to their
!!                         fathers (or the other way around)
!!                         
!! @param[in] fill_thresh  a double precision value containing a threshold for performing amalgamation.
!!                         If the introduced fill-in will be less than
!!                         fill_thresh (in percent) then a node will be
!!                         amalgamated to its father.
!!
subroutine qrm_amalg_tree(n, parent, rowcount, porder, nvar, min_var, fill_thresh, info)

  use qrm_error_mod
  use qrm_mem_mod
  use qrm_common_mod, savesym => qrm_amalg_tree
  implicit none

  integer                      :: n, min_var
  integer                      :: parent(:), rowcount(:), porder(:), nvar(:)
  real(kind(1.d0))             :: fill_thresh
  integer, optional            :: info

  integer, allocatable         :: snodes_ptr(:) ! pointers to the start of a supernode in porder
  integer                      :: nsnodes, i, j, k, f
  integer, allocatable         :: snodes_map(:), snodes_parent(:), kids(:), &
       & snodes_merged(:), snodes_nvar(:), snodes_rc(:)
  real(kind(1.d0))             :: ratio
  integer(kind=8), allocatable :: snodes_fill(:)
  integer(kind=8)              :: totfill, msize, fill
  integer                      :: pv, fpv, rc, frc, col, fcol, tcol, s, ns
  logical                      :: merge

  ! error management
  integer                      :: err
  character(len=*), parameter  :: name='qrm_amalg_tree'

  err = 0
  
  ! sort nodes children in increasing row-count
  call qrm_postorder(parent, n, porder, rowcount, err)
  __QRM_INFO_CHECK(err, name,'qrm_postorder',9999)
  
  ! fundamental supernodes detection
  call qrm_alloc(snodes_ptr, n+1, err)
  call qrm_alloc(kids, n, err)
  __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)
  
  kids = 0
  ! count the number of children for each snode
  do i=1, n
     f = parent(i)
     if(f .ne. 0) kids(f) = kids(f) + 1
  end do

  ! detect snodes
  nsnodes = 1
  snodes_ptr(1) = 1
  do k=2, n
     i = porder(k)
     j = porder(k-1)
     if ( (parent(j) .ne. i) .or. &
        & (rowcount(j) .ne. rowcount(i)+1) .or. &
        & (kids(i) .gt. 1)) then
        nsnodes = nsnodes+1
        snodes_ptr(nsnodes) = k
     end if
  end do
  snodes_ptr(nsnodes+1)=n+1

  ! compress the tree
  ! note that to reduce memory consumption we could
  ! do the amalgamantion on the uncompressed tree TODO
  call qrm_move_alloc(kids, snodes_map)
  call qrm_alloc(snodes_parent, nsnodes, err)
  call qrm_alloc(snodes_nvar, nsnodes, err)
  call qrm_alloc(snodes_rc, nsnodes, err)
  __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)
  
  snodes_nvar = 0

  do i=nsnodes, 1, -1
     pv = porder(snodes_ptr(i))
     snodes_nvar(i) = snodes_ptr(i+1)-snodes_ptr(i)
     snodes_rc(i) = rowcount(pv)
     do j=snodes_ptr(i), snodes_ptr(i+1)-1
        ! map(i) = j means that j belongs to snode i
        snodes_map(porder(j)) = i
     end do
     f = parent(porder(j-1)) ! the father of the last node in the snode
     if(f .eq. 0) then
        snodes_parent(i) = 0
     else
        snodes_parent(i) = snodes_map(f)
     end if
  end do


  allocate(snodes_fill(nsnodes))
  snodes_fill = 0
  call qrm_move_alloc(snodes_map, snodes_merged)
  snodes_merged = 0


  snodes_loop: do s = nsnodes-1, 1, -1

     ! detect the real father of s
     f = snodes_parent(s)
     if(f .ne. 0) then
        do
           if(snodes_merged(f) .eq. 0) exit
           f = snodes_merged(f)
        end do
        
        ! path compression on the chain of amalgamated ancestors
        i = snodes_parent(s)
        do 
           if(i .eq. f) exit
           j = snodes_merged(i)
           snodes_merged(i) = f
           i = j
        end do
     end if

     if(f .ne. s+1) then
        ! don't merge if snode i+1 is not snode i's father
        merge = .false.
        cycle snodes_loop
     end if

     rc   = snodes_rc(s)      ! the row count of the current snode
     frc  = snodes_rc(f)      ! the row count of the father snode

     col  = snodes_nvar(s)    ! the variables in the current snode
     fcol = snodes_nvar(f)    ! the variables in the parent snode
     tcol = col+fcol
     fill = col*(frc+col-rc)  ! the fill added in the case where the nodes are merged
     totfill = snodes_fill(f) + fill

     if(tcol .lt. min_var) then
        merge = .true.
     else if (fill .eq. 0) then
        merge = tcol.lt.min_var*50
     else 
        msize   = (tcol)*(tcol+1)/2 + tcol*(frc-fcol)
        
        ratio   = real(totfill, kind(1.d0))/real(msize,kind(1.d0)) 
        
        merge = (tcol  .lt. min_var*4  .and. ratio .lt. fill_thresh*16) .or.&
             &  (tcol  .lt. min_var*12 .and. ratio .lt. fill_thresh*2) .or. &
             &  (ratio .lt. fill_thresh)

     end if

     if(merge) then
        ! the two nodes ca be merged
        snodes_fill(s)   = totfill
        snodes_nvar(s)   = tcol
        snodes_nvar(f)   = 0
        snodes_merged(f) = s
        snodes_rc(s)     = col+frc
        snodes_rc(f)     = 0
     end if


  end do snodes_loop
  
  nvar=0
  ! rebuild snodes_ptr
  ns = 1
  snodes_ptr(1) = 1
  do i=1, nsnodes
     pv = porder(snodes_ptr(i))
     if(snodes_merged(i) .eq. 0) then
        ns = ns+1
        rowcount(pv) = snodes_rc(i)
        nvar(pv)     = snodes_nvar(i)
        snodes_ptr(ns) = snodes_ptr(ns-1)+snodes_nvar(i)
        if(snodes_ptr(ns).gt.n) exit
     end if
  end do
  ns = ns-1

  call qrm_move_alloc(snodes_merged, snodes_map )
  do i=ns, 1, -1
     pv             = porder(snodes_ptr(i))
     snodes_map(pv) = pv
     f              = parent(porder(snodes_ptr(i+1)-1)) ! the father of the last node in the snode
     do j=snodes_ptr(i)+1, snodes_ptr(i+1)-1
        ! map(i) = j means that j is the principal variable of the snode where i belongs 
        snodes_map(porder(j)) = pv
        rowcount(porder(j))   = 0
        parent(porder(j))     = -pv
     end do
     if(f .eq. 0) then
        parent(pv) = 0
     else
        parent(pv) = snodes_map(f)
     end if
  end do

9999 continue

  call qrm_dealloc(snodes_map)
  call qrm_dealloc(snodes_nvar)
  call qrm_dealloc(snodes_rc)
  call qrm_dealloc(snodes_ptr)
  call qrm_dealloc(snodes_parent)

  if(present(info)) info = err
  return

end subroutine qrm_amalg_tree
