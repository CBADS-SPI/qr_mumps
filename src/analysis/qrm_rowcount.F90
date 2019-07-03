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
!> @file qrm_rowcount.F90
!! This file contains the routine that computes the rowcount for the R factor
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


#include "qrm_common.h"


!> @brief This subroutine computes the rowcount of the R factor 

!> The mothod implemented here is described in:
!!
!! -# J. R. Gilbert, X. S. Li, E. G. Ng, and B. W. Peyton, 2001.
!!     <em>"Computing row and column counts for sparse QR and LU factorization."</em>
!!     BIT Numer. Math. 41, 4, 693--710.
!!
!! -# J. R. Gilbert , E. G. Ng , B. W. Peyton, 1994
!!     <em>"An Efficient Algorithm to Compute Row and Column Counts for Sparse
!!      Cholesky Factorization."</em>
!!     SIAM Journal on Matrix Analysis and Applications,
!!     v.15 n.4, p.1075-1091, Oct. 1994
!!
!! The algorithm on Fig. 3.2 in [1] was (easily) extended with supernodes detection
!! as described in [2].
!!
!! Also, on output, the tree is modified in order to reflect the supernodal structure.
!!
!! @param[in] graph the graph associated to A (or a pruned version if singleton detection was done)
!!
!! @param[out] rc the row count. rc(i)=k means that in the k-th row of
!!                R there are k nonzeroes; k=0 for all the
!!                subordinate variables.  This also gives us the
!!                size of the rows of front i.
!!
!! @param[in,out] porder on input an integer array containing a
!!                       postorder of the tree. On output it contains
!!                       an equivalent postorder where principal
!!                       variables come always before the correspondig
!!                       subordinates. Other routines rely on this
!!                       postorder and thus it should never be
!!                       changed.
!!
!! @param[in,out] parent an integer array containing the elimination tree in input
!!                       and the assembly tree on output. The meaning of parent on output is:
!!                       -# parent(i) = j>0: i is the principal
!!                                           variable of a node and j
!!                                           if the principal variable
!!                                           of its father node
!!                       -# parent(i) = j=0: i is a principal variable of a root node
!!                       -# parent(i) = j<0: i is a subordinate
!!                                           variable inside a node
!!                                           whose principal variable
!!                                           is j.
!!
!! Example output:
!! @verbatim
!!         +---+
!!         |7  |
!!         | 6 |           
!!         |  5|           parent=(/ -1, 7, -4, 7, -7, -7, 0 /)
!!         +---+
!!        /     \
!!       /       \
!!   +--+         +--+
!!   |2 |         |4 |
!!   | 1|         | 3|
!!   +--+         +--+
!!@endverbatim
!!
!! @warning at the moment the detection of supernodes is turned off (see comments in the code).
!!          This is due to the fact that supernodes are detected by the amalgamation routine
!!          @link qrm_amalg_tree @endlink. As a result, all the entries of parent will be
!!          positive on output.

subroutine _qrm_rowcount(graph, parent, porder, rc, info)


  use _qrm_spmat_mod
  use qrm_mem_mod
  use qrm_common_mod
  implicit none

  type(_qrm_spmat_type)           :: graph
  integer                         :: parent(:), porder(:), rc(:)
  integer, optional               :: info
  
  ! level    : it is the level of node i, i.e., the distance from its root
  ! fst_desc : it is the first descendent of node i according to the postorder

  integer, allocatable            :: fst_desc(:), stack(:)
  integer, allocatable            :: first(:), hptr(:), hjcn(:), rcnt(:), last(:), ipord(:)
  integer, allocatable            :: prev_f(:), prev_nbr(:), setpath(:)
  integer                         :: i, curr, fd, hp, rlev, j
  integer                         :: f, ptr, k, ii, u, ref, p_leaf, q, jj

  ! error management
  integer                         :: err
  character(len=*), parameter     :: name='qrm_rowcount'

  info = 0
  
  call qrm_alloc(fst_desc, graph%n, err)
  __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)
  call qrm_alloc(stack   , graph%n, err)
  __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)
  
  rc = 0
  fst_desc=-1
  do i=1, graph%n
     curr = porder(i)
     fd   = curr
     if(fst_desc(curr) .eq. -1) then
        ! this is a leaf node. Initialize its row_count to 1
        ! see [1], Fig 3.2
        rc(curr) = 1
     end if

     ! now we start going up thre tree seeting up fd as first
     ! descendent of all its ancestors
     do
        if (fst_desc(curr) .gt. 0) exit
        fst_desc(curr) = fd
        if (parent(curr) .eq. 0) exit
        curr = parent(curr)
     end do
  end do

  ! now we need to compute the "higher adjacency sets" (see [1] pag. 704)
  ! Assume first(i)=j means that the first nonzero in row i has column index j.
  ! The higher adjaceny set of column j is defined as the union of the column indices
  ! of all the rows i for which first(i)=j.
  ! The computation of the higher adjacency sets may be expensive in memory or in perf.
  ! Here we choose (for the moment) to spend some more time on it but have a better
  ! memory efficiency.
  call qrm_move_alloc(stack, last)
  if(err.eq.0) call qrm_alloc(ipord, graph%n  , err)
  if(err.eq.0) call qrm_alloc(hptr , graph%n+1, err)
  if(err.eq.0) call qrm_alloc(rcnt , graph%n  , err)
  if(err.eq.0) call qrm_alloc(first, graph%m  , err)
  __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)

  ! compute inverse postorder
  do j=1, graph%n
     ipord(porder(j)) = j
  end do

  first = 0
  rcnt  = 0
  last  = 0
  ! we make a first pass to count the number of elements in each higher adjacency set
  do k=1, graph%n
     j = porder(k)
     do ii=graph%jptr(j), graph%jptr(j+1)-1
        i = graph%irn(ii)
        f = first(i)
        if(f .eq. 0) then
           first(i) = j
        else if (k .gt. last(f)) then
           rcnt(f) = rcnt(f)+1
           last(f) = k
        end if
     end do
  end do

  ! transform the counts into pointers
  hptr(1) = 1
  do k=1, graph%n
     hptr(k+1) = hptr(k)+rcnt(k)
  end do

  ! allocate hjcn
  call qrm_alloc(hjcn, hptr(graph%n+1), err)
  __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)

  rcnt  = 0
  last  = 0
  ! second loop to fill
  do k=1, graph%n
     j = porder(k)
     do ii=graph%jptr(j), graph%jptr(j+1)-1
        i = graph%irn(ii)
        f = first(i)
        if ( (k .gt. ipord(f)) .and. (k .gt. last(f)) ) then
           ptr = hptr(f) + rcnt(f)
           hjcn(ptr) = j
           rcnt(f) = rcnt(f)+1
           last(f) = k
        end if
     end do
  end do


  ! at this point we have everything we need to compute the rowcount.
  ! This is the algorithm on Figure 3.2 in [1]
  call qrm_move_alloc(last, prev_f)
  call qrm_move_alloc(rcnt, prev_nbr)
  call qrm_dealloc(first)
  ! setpath is used for doing path compression
  call qrm_alloc(setpath, graph%n, err)
  __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)

  do j=1, graph%n
     setpath(j)=j
  end do
  prev_f   = 0
  prev_nbr = 0
  
  ! during the symbolic facto we flag as negative all the first
  ! variables of a supernode. At the end, we make a simple loop
  ! to modify the tree so that it reflects the supernodal structure
  
  do jj=1, graph%n
     j = porder(jj)
     if(parent(j) .ne. 0) then
        rc(parent(j)) = rc(parent(j))-1
     else
        ! if porder(jj) is a root, porder(jj+1) is the beginning of a supernode
        if(jj .ne. graph%n) porder(jj+1) = porder(jj+1)
     end if
     
     do k=hptr(j), hptr(j+1)-1
        u = hjcn(k)
        if(prev_nbr(u) .eq. 0) then
           ref = 0
        else
           ref = ipord(prev_nbr(u))
        end if

        if(ipord(fst_desc(j)) .gt. ref) then
           rc(j)  = rc(j)+1
           p_leaf = prev_f(u)
           ! porder(jj) is a leaf or a row subtree, thus the first variable in a
           ! supernode. flag it
           if (p_leaf .ne. 0) then
              q = setfind(setpath, p_leaf)
              ! q is a leaf or a row subtree, thus the first variable in a
              ! supernode. flag it
              rc(q) = rc(q) -1
           end if
           prev_f(u) = j
        end if
        prev_nbr(u) = j
     end do
     call setunion(setpath, j, parent(j))

  end do

  do jj=1, graph%n-1
     j=porder(jj)
     if(parent(j) .ne. 0) rc(parent(j)) = rc(parent(j)) + rc(j)
  end do

  ! WARNING: this commented code is meant to detect supernodes in the tree based
  ! on how the nodes where flagged above. It is currently replaced by the small
  ! loop below as supernodes are detected in the amalgamation routine. As a result,
  ! on output all the entries of parent will be positive.

  ! ! now we have to update the tree in order to reflect the supernodal structure.
  ! ! we simply traverse the postorder looking for the variables we flagged above
  ! i=1
  ! do while (i .lt. graph%n)
     ! porder(i) = unflip(porder(i))
     ! j = i+1
     ! do while( (porder(j) .gt. 0) .and. (j .le. graph%n))
        ! j = j+1
     ! end do
     ! ! at this point j is the first variable of the next supernode and
     ! ! j-1 is the principal variable of the supernode starting at i.
     ! ! All the variables i:j-2 are subordinate to j-1
     ! parent(porder(i)) = parent(porder(j-1))
     ! do k=i+1, j-1
        ! parent(porder(k)) = -porder(i)
        ! rc(porder(k)) = 0
     ! end do
     ! ! we want the supervariable to be in first position in the porder
     ! i = j
  ! end do

  ! unflip the last in case it's a supervariable
  ! porder(graph%n) = unflip(porder(graph%n))

  do i=1, graph%n
     f = parent(i)
     if (f .gt. 0) then
        if(parent(f) .lt. 0) parent(i) = -parent(f)
     end if
  end do

9999 continue
  ! cleanup and return
  call qrm_dealloc(setpath )
  call qrm_dealloc(hptr    )
  call qrm_dealloc(hjcn    )
  call qrm_dealloc(prev_nbr)
  call qrm_dealloc(prev_f  )
  call qrm_dealloc(fst_desc)
  call qrm_dealloc(ipord   )

  if(present(info)) info = err
  return

contains

  function setfind(setpath, p_leaf)
    implicit none
    integer :: setpath(:), p_leaf, setfind

    integer :: q, c, tmp

    q=p_leaf

    do while (setpath(q) .ne.q)
       q = setpath(q)
    end do

    ! path compression
    c = p_leaf
    do while (c .ne.q)
       tmp = setpath(c)
       setpath(c) = q
       c = tmp
    end do

    setfind = q
    return
  end function setfind


  subroutine setunion(setpath, j, pj)
    implicit none
    integer :: setpath(:), j, pj
    if(pj .ne. 0) setpath(j) = pj
    return
  end subroutine setunion


  function flip(p)
    integer :: p, flip
    if(p .gt. 0) then
       flip = -p
    else
       flip = p
    end if
  end function flip
    
  function unflip(p)
    integer :: p, unflip
    if(p .lt. 0) then
       unflip = -p
    else
       unflip = p
    end if
  end function unflip

end subroutine _qrm_rowcount
