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
!> @file qrm_ata_graph.F90
!! This file contains the routine that computes the graph of A^T * A
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief This subroutine computes the fill reducing ordering using METIS
!! 
!! @param[in] g_csc  the graph (in CSC format) associated to A
!! 
!! @param[out] ata_graph  the graph of A'*A in CSR format
!! 
subroutine _qrm_ata_graph(g_csc, ata_graph, info)
  use _qrm_spmat_mod
  use qrm_error_mod
  use qrm_mem_mod

  implicit none

  type(_qrm_spmat_type), intent(in)  :: g_csc
  type(_qrm_spmat_type), intent(out) :: ata_graph
  integer, optional                  :: info
  
  type(_qrm_spmat_type)              :: g_csr
  integer                            :: i, j, row1, col1, row2, col2
  integer, allocatable               :: mark(:)

  ! error management
  integer                            :: err, err2
  character(len=*), parameter        :: name='qrm_ata_graph'
  
  err = 0; err2 = 0
  
  call _qrm_spmat_convert(g_csc, g_csr, 'csr', values=.false., info=err)
  __QRM_INFO_CHECK(err, name, 'qrm_spmat_convert', 9999)

  call qrm_alloc(ata_graph%iptr,g_csc%n+2, err)
  __QRM_INFO_CHECK(err, name, 'qrm_alloc', 9999)
  ata_graph%iptr = 0
  ata_graph%iptr(1:2) = 1

  call qrm_alloc(mark, g_csc%n, err)
  __QRM_INFO_CHECK(err, name, 'qrm_alloc', 9999)

  mark = 0

  do col1=1, g_csc%n
     ! loop over all the columns of A
     do i=g_csc%jptr(col1), g_csc%jptr(col1+1)-1
        ! for each nnz in the column, we go through the corresponding
        ! row and count all the i,j pairs
        row1 = g_csc%irn(i)

        do j=g_csr%iptr(row1), g_csr%iptr(row1+1)-1
           col2 = g_csr%jcn(j)
           ! now, the element (col1,col2) is present in A'*A. Check if
           ! it was already found, otherwise add it.
           ! skip the diagonal
           if(col1 .eq. col2) cycle
           if(mark(col2) .lt. col1) then
              mark(col2) = col1
              ata_graph%iptr(col1+2) = ata_graph%iptr(col1+2)+1
           end if
        end do
     end do
  end do

  do i=3, g_csc%n+2
     ata_graph%iptr(i) = ata_graph%iptr(i)+ata_graph%iptr(i-1)
  end do

  ata_graph%nz = ata_graph%iptr(g_csc%n+2)
  call qrm_alloc(ata_graph%jcn, ata_graph%nz, err)
  __QRM_INFO_CHECK(err, name, 'qrm_alloc', 9999)

  ! second pass to fill up
  mark = 0
  do col1=1, g_csc%n
     ! loop over all the columns of A
     do i=g_csc%jptr(col1), g_csc%jptr(col1+1)-1
        ! for each nnz in the column, we go through the corresponding
        ! row and count all the i,j pairs
        row1 = g_csc%irn(i)

        do j=g_csr%iptr(row1), g_csr%iptr(row1+1)-1
           col2 = g_csr%jcn(j)
           ! now, the element (col1,col2) is present in A'*A. Check if
           ! it was already found, otherwise add it.
           ! skip the diagonal
           if(col1 .eq. col2) cycle
           if(mark(col2) .lt. col1) then
              mark(col2) = col1
              ata_graph%jcn(ata_graph%iptr(col1+1)) = col2
              ata_graph%iptr(col1+1) = ata_graph%iptr(col1+1)+1
           end if
        end do
     end do
  end do

  ata_graph%n = g_csc%n
  ata_graph%m = g_csc%n

  call _qrm_spmat_destroy(g_csr, all=.true., info=err)
  __QRM_INFO_CHECK(err, name, 'qrm_spmat_destroy', 9999)
  call qrm_dealloc(mark, err)
  __QRM_INFO_CHECK(err, name, 'qrm_dealloc', 9999)

  if(present(info)) info = err
  return
  
9999 continue
  ! cleanup and return
  call _qrm_spmat_destroy(g_csr, all=.true.)
  call _qrm_spmat_destroy(ata_graph, all=.true.)
  call qrm_dealloc(mark)

  if(present(info)) info = err
  return

end subroutine _qrm_ata_graph
