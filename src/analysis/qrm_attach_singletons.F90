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
!> @file qrm_attach_singletons.F90
!! This file contains a routine that merges the results of the singletons detection into the
!! results of the analysis phase
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief This subroutine merges the results of the singletons detection into the
!! results of the analysis phase

!> Once the singletons are detected in the @link
!! _qrm_detect_singletons @endlink routine, the matrix graph is
!! purged, i.e., the rows and columns corresponding to the singletons
!! are removed. The subsequent operations (column ordering,
!! elimination tree, symboli factorization etc.) are executed on the
!! purged graph. Thus, at some point, the results of all these
!! operations have to be merged with what was found at the singleton
!! detection. Specifically all the singeltons are grouped into a front
!! (number 1) which will be recognized as its row count is set to
!! -1. All the adata memebers are modified accordingly in this routine.
!!
!! @param[in,out] qrm_mat The instance of the problem where the
!!                        singletons have to be attached. Only the
!!                        adata member will be modified (in all its
!!                        members)
!! 
!! @param[in] scol the list of column singletons
!! 
!! @param[in] srow the list of row singletons
!! 
!! @param[in] mrperm an integer array storing the mapping to the row
!!                    indices of the purged matrix. mrperm(j) = k
!!                    means that column j of the purged graph
!!                    corresponds to column k in the original matrix.
!!                    A row permutation is later computed on the
!!                    purged matrix and, thus, mrperm will be used to
!!                    map this permutation back on the original
!!                    matrix.
!!
!! @param[in] mcperm an integer array storing the mapping to the row
!!                    indices of the purged matrix. mcperm(j) = k
!!                    means that row j of the purged graph corresponds
!!                    to row k in the original matrix.  Because the
!!                    ordering method will be executed on the purged
!!                    graph, this array is necessary to map the
!!                    computed ordering onto the columns of A.
!! 
!! @param[in] nrsing the number of row singletons
!! 
!! @param[in] ncsing the number of column singletons
!! 
subroutine _qrm_attach_singletons(qrm_mat, scol, srow, mrperm, mcperm, nrsing, ncsing)
  
  use _qrm_spmat_mod
  use qrm_adata_mod
  implicit none

  type(_qrm_spmat_type) :: qrm_mat
  integer                :: scol(:), srow(:), mcperm(:), mrperm(:)
  integer                :: nrsing, ncsing


  integer, allocatable :: tmp(:)
  integer :: m, n, i, nnodes

  ! restore the values of m and n (to be equal to the original matrix
  ! dimensions)
  m = qrm_mat%m+nrsing 
  n = qrm_mat%n+ncsing

  ! restore the values of nnodes. One extra node will be attached to
  ! the tree corresponding to the set of singletons
  nnodes = qrm_mat%adata%nnodes+1


  ! attach the singleton columns to the column permutation
  call qrm_alloc(tmp, n)
  tmp(1:ncsing) = scol
  do i=1, qrm_mat%n
     tmp(ncsing+i) = mcperm(qrm_mat%adata%cperm(i))
  end do
  call qrm_dealloc(qrm_mat%adata%cperm)
  call move_alloc(tmp, qrm_mat%adata%cperm)



  ! attach the singleton rows to the row permutation
  call qrm_alloc(tmp, m)
  tmp(1:nrsing) = srow
  do i=1, qrm_mat%m
     tmp(nrsing+i) = mrperm(qrm_mat%adata%rperm(i))
  end do
  call qrm_dealloc(qrm_mat%adata%rperm)
  call move_alloc(tmp, qrm_mat%adata%rperm)



  ! attach a new node to the cp_ptr array
  call qrm_alloc(tmp, nnodes+1)
  tmp(1)=1
  tmp(2:nnodes+1) = qrm_mat%adata%cp_ptr(1:qrm_mat%adata%nnodes+1)+ncsing
  call qrm_dealloc(qrm_mat%adata%cp_ptr)
  call move_alloc(tmp, qrm_mat%adata%cp_ptr)



  ! attach a new node to the parent array
  call qrm_alloc(tmp, nnodes)
  tmp(1)=0
  do i=1, qrm_mat%adata%nnodes
     if(qrm_mat%adata%parent(i) .eq. 0) then
        tmp(i+1)=0
     else
        tmp(i+1) = qrm_mat%adata%parent(i)+1
     end if
  end do
  call qrm_dealloc(qrm_mat%adata%parent)
  call move_alloc(tmp, qrm_mat%adata%parent)



  ! attach a new node to the rc array
  call qrm_alloc(tmp, nnodes)
  tmp(1)=-1
  tmp(2:nnodes) = qrm_mat%adata%rc(1:qrm_mat%adata%nnodes)
  call qrm_dealloc(qrm_mat%adata%rc)
  call move_alloc(tmp, qrm_mat%adata%rc)



  ! attach a new node to the stair array
  call qrm_alloc(tmp, nnodes)
  tmp(1)=nrsing
  tmp(2:nnodes) = qrm_mat%adata%stair(1:qrm_mat%adata%nnodes)+nrsing
  call qrm_dealloc(qrm_mat%adata%stair)
  call move_alloc(tmp, qrm_mat%adata%stair)



  ! attach a new node to the childptr array and fix the child array
  call qrm_alloc(tmp, nnodes+1)
  tmp(1)=1
  tmp(2:nnodes+1) = qrm_mat%adata%childptr(1:qrm_mat%adata%nnodes+1)
  qrm_mat%adata%child = qrm_mat%adata%child+1
  call qrm_dealloc(qrm_mat%adata%childptr)
  call move_alloc(tmp, qrm_mat%adata%childptr)



  ! attach a new node to the nfrows array
  call qrm_alloc(tmp, nnodes)
  tmp(1)=nrsing
  tmp(2:nnodes) = qrm_mat%adata%nfrows(1:qrm_mat%adata%nnodes)
  call qrm_dealloc(qrm_mat%adata%nfrows)
  call move_alloc(tmp, qrm_mat%adata%nfrows)



  ! attach a new node to the fcol and fcol_ptr arrays
  call qrm_alloc(tmp, nnodes+1)
  tmp(1)=1
  tmp(2:nnodes+1) = qrm_mat%adata%fcol_ptr(1:qrm_mat%adata%nnodes+1)
  call qrm_dealloc(qrm_mat%adata%fcol_ptr)
  call move_alloc(tmp, qrm_mat%adata%fcol_ptr)
  do i=1, qrm_mat%adata%fcol_ptr(nnodes+1)-1
     qrm_mat%adata%fcol(i) = mcperm(qrm_mat%adata%fcol(i))
  end do



  ! fix the leaves and small arrays
  qrm_mat%adata%leaves = qrm_mat%adata%leaves+1
  call qrm_alloc(tmp, nnodes)
  tmp(1)=0
  tmp(2:nnodes) = qrm_mat%adata%small(1:qrm_mat%adata%nnodes)
  call qrm_dealloc(qrm_mat%adata%small)
  call move_alloc(tmp, qrm_mat%adata%small)
  
  qrm_mat%adata%nnodes = nnodes
  qrm_mat%adata%ncsing = ncsing
  qrm_mat%adata%nrsing = nrsing


  return


end subroutine _qrm_attach_singletons
