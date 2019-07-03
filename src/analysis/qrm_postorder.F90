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
!> @file qrm_postorder.F90
!! This file contains the routine that computes a postorder traversal of a tree
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!>   @brief This subroutine computes a postorder by traversing a tree in dfs.
!!   
!!   @param[in] parent  integer array of size n. parent(i)=j means that node j is the
!!                      father of node i in the tree
!!   
!!   @param[in] n       number of nodes in the tree
!!   
!!   @param[out] porder an integer array of size n containing the postorder
!!
!!   @param[in] weight  an optional array containing nodes weights. If present, the
!!                      children of each node will be sorted by increasing weight.
!!   
subroutine qrm_postorder(parent, n, porder, weight, info)

  use qrm_mem_mod
  implicit none

  integer                         :: n
  integer                         :: parent(:), porder(:)
  integer, optional               :: weight(:)
  integer, optional               :: info
  
  integer, allocatable            :: son(:), brother(:), stack(:)
  integer                         :: i, father, br, head, hp, pp, t, w, next

  ! error management
  integer                         :: err
  character(len=*), parameter     :: name='qrm_postorder'

  err = 0

  if(err.eq.0) call qrm_alloc(son,     n, err)
  if(err.eq.0) call qrm_alloc(brother, n, err)
  if(err.eq.0) call qrm_alloc(stack,   n, err)
  __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)

  son = 0

  ! build a tree that can be traversed top-to-bottom
  if(present(weight)) then
     ! use stack as a workspace
     stack = 0
     do i=1, n
        w = weight(i)
        brother(i) = stack(w)
        stack(w)   = i
     end do

     do w=n, 1, -1
        i = stack(w)
        do
           if(i .eq. 0) exit
           next   = brother(i)
           father = parent(i)
           if(father .ne. 0) then
              brother(i) = son(father)
              son(father) = i
           end if
           i = next
        end do
     end do
  else
     do i=n, 1, -1
        father = parent(i)
        if (father .ne. 0) then
           br          = son(father)
           brother(i)  = br
           son(father) = i
        end if
     end do
  end if


  head = 0
  hp   = 0
  pp   = 1
  ! the tree is processed in dfs. Starting from a root, we go down
  ! and put all the encountered nodes on a stack. When we reach the bottom,
  ! we pop the leaf from the stack, we put it in the postorder and we go down again
  ! along the branch starting from the brother of the node just popped
  do t=1, n
     if (parent(t) .ne. 0) cycle
     ! at this point t is the root of a tree
     hp        = hp+1
     stack(hp) = t
     head      = t
     do
        if(son(head) .eq. 0) then
           ! we reached the bottom
           porder(pp) = head
           pp = pp+1
           hp = hp-1
           if (parent(head) .ne. 0) then
              son(parent(head)) = brother(head)
           end if
           if (hp .eq. 0) then
              ! tree is exhausted
              exit
           end if
           head = stack(hp)
        else
           ! go down one more level
           hp = hp+1
           stack(hp) = son(head)
           head = son(head)
        end if
     end do
  end do

  
9999 continue
  ! cleanup and return
  call qrm_dealloc(son)
  call qrm_dealloc(brother)
  call qrm_dealloc(stack)

  if(present(info)) info = err
  return

end subroutine qrm_postorder
