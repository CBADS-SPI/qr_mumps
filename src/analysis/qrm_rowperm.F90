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
!> @file qrm_rowperm.F90
!! this file contains the routine which computes a row permutation that puts the input
!! matrix into a staircase format.
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief This routine computes a row permutation that puts the input
!! matrix into a
!! staircase format.

!> This subroutine computes a row permutation rperm such that
!! A(rperm, cperm) is as much as possible upper triangular.
!! The idea is to have someting like:
!! @verbatim
!!   |x       |
!!   |x       |     in this case stair would be:
!!   | x      |     stair=(/ 2, 5, 7, 7, 9 /)
!!   | x      |
!!   | x      |
!!   |  xx    |
!!   |  xx    |
!!   |    x   |
!!   |    x   |
!! @endverbatim
!!
!! @param[in] graph  a qrm_spmat_type containing the matrix
!! 
!! @param[in] cperm  a column permutation
!! 
!! @param[out] rperm the computed row permutation
!!
!! @param[in] nvar   nvar(i) is the number of variables in the node whose principal variable is i
!! 
!! @param[out] stair stair(i) is the number of rows in step i of the
!!                   stair structure computed on the input graph
!!
subroutine _qrm_rowperm(graph, cperm, rperm, nvar, stair, info)

  use _qrm_spmat_mod
  use qrm_mem_mod
  use qrm_common_mod
  implicit none

  type(_qrm_spmat_type)           :: graph
  integer                         :: cperm(:), rperm(:), stair(:), nvar(:)
  integer, optional               :: info
  
  integer, allocatable            :: mark(:)
  integer                         :: ii, jj, i, j, node, nv, pnt, rpnt

  ! error management
  integer                         :: err
  character(len=*), parameter     :: name='qrm_rowperm'
  
  err = 0
  
  call qrm_alloc(mark, graph%m, err)
  __QRM_INFO_CHECK(err, name, 'qrm_alloc', 9999)
  
  mark  = 0
  stair = 0
  ! scan the matrix by columns

  rpnt = 0
  pnt=1
  do 
     if (pnt .gt. graph%n) exit
     node = cperm(pnt)
     stair(node) = rpnt
     nv = nvar(node)
     do jj = pnt, pnt+nv-1
        j = cperm(jj)
        do ii=graph%jptr(j), graph%jptr(j+1)-1
           i = graph%irn(ii)
           if(mark(i) .eq. 0) then
              stair(node) = stair(node)+1
              rpnt = rpnt+1
              rperm(rpnt) = i
              mark(i) = j
           end if
        end do
     end do
     pnt = pnt+nv
     if (pnt .gt. graph%n) exit
  end do

  ! look for empty rows and put them at the end
  do i=1, graph%m
     if(mark(i) .eq. 0) then
        rpnt = rpnt+1
        rperm(rpnt) = i
        mark(i) = i
     end if
  end do

  call qrm_dealloc(mark, err)
  __QRM_INFO_CHECK(err, name, 'qrm_dealloc', 9999)

9999 continue
  if(present(info)) info = err
  return

end subroutine _qrm_rowperm
