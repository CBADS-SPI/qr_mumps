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
!> @file qrm_elim_tree.F90
!! This file contains the routine that computes the eliimation tree
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


#include "qrm_common.h" 


!> @brief This subroutine builds the elimination
!! tree for A'A. 

!> The graph of A is contained in the "graph" input argument and A's
!> columns are permuted accoridng cperm. This is an implementation of
!> the algorithm described in
!!
!! Gilbert, J. R., Li, X. S., Ng, E. G., and Peyton, B. W. 2001.
!! <em>"Computing row and column counts for sparse QR and LU factorization."</em>
!! BIT Numer. Math. 41, 4, 693--710.
!!
!! @param[in] graph    the graph associated to matrix A in CSC format
!!
!! @param[in] cperm    the pivotal order (i.e., the column permutation previously
!!                     computed by the @link _qrm_do_ordering @endlink routine
!! 
!! @param[out] parent  the elimitation tree. parent(i)=j means that
!!                     node j is the parent of node i in the
!!                     elimination tree. The parent of a root node is
!!                     equal to 0. This tree is of size n where n is
!!                     the size of the matrix
!!
subroutine _qrm_elim_tree(graph, cperm, parent, info)

  use _qrm_spmat_mod
  use qrm_mem_mod
  implicit none

  type(_qrm_spmat_type), intent(in)  :: graph
  integer, intent(in)                :: cperm(:)
  integer, allocatable               :: parent(:)
  integer, optional                  :: info

  integer, allocatable               :: ancestor(:), prev_col(:)
  integer                            :: iidx, jidx, i, j, k

  ! error management
  integer                            :: err, err2
  character(len=*), parameter        :: name='qrm_elim_tree'

  err = 0

  call qrm_alloc(ancestor, graph%n, err)
  __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)
  call qrm_alloc(prev_col, graph%m, err)
  __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)
  call qrm_realloc(parent, graph%n, err) ! realloc parent in case it was not allocated before
  __QRM_INFO_CHECK(err, name,'qrm_realloc',9999)

  ! ancestor is used for doing path compression
  ancestor=0
  prev_col=0
  parent  =0

  do jidx = 1, graph%n
     j = cperm(jidx)
     do iidx = graph%jptr(j), graph%jptr(j+1)-1
        i = graph%irn(iidx)
        k = prev_col(i)
        if(k .ne. 0) call add_node(ancestor, parent, k, j)
        prev_col(i) = j
     end do
  end do

9999 continue
  ! cleanup and return
  err2 = 0
  call qrm_dealloc(ancestor, err2)
  if(err2.eq.0) call qrm_dealloc(prev_col, err2)

  if(present(info)) info = merge(err, err2, err.ne.0)
  return

contains
  
  subroutine add_node(ancestor, parent, k, j)
    implicit none
    integer, allocatable               :: parent(:)
    integer, allocatable               :: ancestor(:)
    integer                            :: k, j

    integer                            :: r

    continue
    do
       r = ancestor(k)

       ancestor(k)=j

       if (r .eq. j) then
          ! node was already added, do nothing
          return
       end if

       ancestor(k)=j

       if(r .eq. 0) then
          ! we reached the top of the subtree
          ! add j on top
          parent(k) = j
          return
       end if

       k = r

    end do


  end subroutine add_node

end subroutine _qrm_elim_tree
