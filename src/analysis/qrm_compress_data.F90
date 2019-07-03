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
!> @file qrm_compress_data.F90
!! this file contains the subroutine that compresses the result of several operations done
!! during the analysis to a size that is proportional to the number of nodes in the elimination
!! tree. These data are of size ~n in input
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


#include "qrm_common.h"


!> @brief This routine compresses the results of a number of operations in the analysis phase.
!! Basically, the input data is of size n and the output of size adata%nnodes (which
!! is the number of nodes in the elimination tree).

!! The arrays on input are all of size n (where n is the column size of
!! the original matrix). The purpose of this subroutine is to store the
!! same information into arrays of size nsteps, where nsteps is the
!! number of nodes in the assembly tree. Only porder will remain of size n
!!
!! @param[out] adata  on output the following fields will be
!!          modified:
!!          - cp_ptr
!!          - rc
!!          - parent
!!          - nnodes
!!          - stair
!!          - child. for each node, it contains the list of its children
!!          - childptr. pointers to the list of children. chil(childptr(i),...,childptr(i+1)-1)
!!            contains all the children of node i
!!          - icperm. the inverse column permutation is built
!!          
!! @param[in] porder  contains the postorder. porder(i)=k means that the i-th column
!!                    in the computed ordering is column k in the original matrix.
!!                    Inthis postorder principal variables always come before the
!!                    corresponding subordinate variables.
!!          
!! @param[in] parent  contains the assembly tree. parent(i)=k:
!!                    - k>0: means that i is a principal variable in a node and k
!!                      is the principal variable in the father's node
!!                    - k<0: means that i is a subordinate variable in a supernode
!!                      whose principal variable is k.
!!          
!! @param[in] rc      this array contains the rowcount, i.e., rc(i)=k means that
!!                    in the R factor the rows corresponding to the node whose
!!                    principal variable is i have k nonzeroes. rc(i)=0 for all
!!                    subordinate varibales and rc(i)=-1 for all the column singletons
!!          
!! @param[in] stair   stair(i) contains the number of rows in the step related to
!!                    node i. see @link ::_qrm_rowperm_ @endlink
!!          
!! @param[in] n       the number of columns in the original matrix (i.e. the size of the
!!                    input arrays)
!!          
subroutine qrm_compress_data(adata, porder, parent, rc, stair, n, info)
  
  use qrm_mem_mod
  use qrm_adata_mod
  use qrm_common_mod, protect => qrm_compress_data
  implicit none
  
  type(qrm_adata_type)            :: adata
  integer                         :: porder(:), parent(:), rc(:), stair(:)
  integer                         :: n
  integer, optional               :: info
  
  integer                         :: i, pnt, svar, f, ss, cnt
  integer, allocatable            :: work(:)
  
  ! error management
  integer                         :: err, err2
  character(len=*), parameter     :: name='qrm_compress_data'

  err = 0; err2 = 0

  adata%nnodes = 0
  do i=1, n
     if(rc(i) .ne. 0) adata%nnodes = adata%nnodes+1
  end do

  call qrm_alloc(adata%cp_ptr, adata%nnodes+1+1, err)
  __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)
  call qrm_alloc(work, n, err)
  __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)
  
  work = 0
  adata%cp_ptr=0
  ! build pointers and inverse mapping
  ! work(i)=k means that principal variable i is in node number k
  adata%cp_ptr(1) = 1
  work(porder(1)) = 1
  pnt             = 2
  do i=2, n
     svar=porder(i)
     if (rc(svar) .ne. 0) then
        adata%cp_ptr(pnt) = i
        work(svar) = pnt
        pnt = pnt+1
     end if
  end do
  adata%cp_ptr(adata%nnodes+1) = n+1

  ! build adata%parent
  ! build adata%rc
  ! build adata%child and adata%childptr
  if(err.eq.0) call qrm_alloc(adata%parent,   adata%nnodes+1,   err)
  if(err.eq.0) call qrm_alloc(adata%rc,       adata%nnodes+1,   err)
  if(err.eq.0) call qrm_alloc(adata%stair,    adata%nnodes+1,   err)
  if(err.eq.0) call qrm_alloc(adata%child,    adata%nnodes+1,   err)
  if(err.eq.0) call qrm_alloc(adata%childptr, adata%nnodes+1+1, err)
  __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)

  adata%childptr = 0
  ss = 0
  do i=1, adata%nnodes
     svar = porder(adata%cp_ptr(i))
     f    = parent(svar)
     adata%rc(i) = rc(svar)
     ! adata%nfrows(i) = stair(svar)-ss
     adata%stair(i) = stair(svar)
     ss = adata%stair(i)
     if(f .eq. 0) then
        adata%parent(i) = 0
     else     
        adata%parent(i) = work(f)
        adata%childptr(work(f)+1) = adata%childptr(work(f)+1)+1
     end if
  end do

  adata%childptr(1) = 1
  do i=2, adata%nnodes+1
     adata%childptr(i) = adata%childptr(i-1)+adata%childptr(i)
  end do


  adata%child=0
  work(1:adata%nnodes)=0
  do i=1, adata%nnodes
     f = adata%parent(i)
     if (f .ne. 0) then
        adata%child(adata%childptr(f)+work(f)) = i
        work(f) = work(f)+1
     end if
  end do

  ! add a fake node, father of all the roots. This guarantees that the
  ! tree is always a tree and not a forest
  adata%cp_ptr(adata%nnodes+2) = n+1
  cnt = adata%childptr(adata%nnodes+1)
  do i=1, adata%nnodes
     if(adata%parent(i).eq.0) then
        adata%parent(i) = adata%nnodes+1
        adata%child(cnt) = i
        cnt = cnt+1
        adata%childptr(adata%nnodes+2) = cnt
     end if
  end do
  adata%parent(adata%nnodes+1) = 0
  adata%rc(adata%nnodes+1) = 0
  adata%stair(adata%nnodes+1) = adata%stair(adata%nnodes)

  adata%nnodes = adata%nnodes+1
  
  call qrm_dealloc(work, err)
  __QRM_INFO_CHECK(err, name,'qrm_dealloc',9999)

  if(present(info)) info=err
  return

9999 continue
  ! cleanup and return
  if(err2.eq.0) call qrm_dealloc(work,           err2)
  if(err2.eq.0) call qrm_dealloc(adata%parent,   err2)
  if(err2.eq.0) call qrm_dealloc(adata%rc,       err2)
  if(err2.eq.0) call qrm_dealloc(adata%stair,    err2)
  if(err2.eq.0) call qrm_dealloc(adata%child,    err2)
  if(err2.eq.0) call qrm_dealloc(adata%childptr, err2)

  if(present(info)) info=merge(err, err2, err.ne.0)
  return

end subroutine qrm_compress_data
