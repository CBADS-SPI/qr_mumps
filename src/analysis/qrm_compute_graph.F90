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
!> @file qrm_compute_graph.F90
!! This file contains the routine that computes the graph of a matrix
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief Computes the adjacency graph of a matrix

!> This subroutine computes the column graph associated to an input matrix qrm_mat 
!! COO format. The output graph has no duplicates as well as no self-edges
!!
!! @param[in] qrm_mat the input matrix
!!
!! @param[out] graph  the adjacency graph in CSC format
!!
subroutine _qrm_compute_graph(qrm_mat, graph, transp, info)

  use _qrm_spmat_mod
  use qrm_error_mod
  use qrm_mem_mod
  use _qrm_analysis_mod, savesym => _qrm_compute_graph

  implicit none

  type(_qrm_spmat_type)              :: qrm_mat
  type(_qrm_spmat_type), intent(out) :: graph
  character                          :: transp
  integer, optional                  :: info

  integer                            :: i, j, pnt, savepnt, dups, ii
  integer, allocatable               :: dupsmap(:)

  ! error management
  integer                            :: err, err2
  character(len=*), parameter        :: name='qrm_compute_graph'

  err = 0; err2 = 0
  
  ! At this moment, if the matrix is centralized, we just need to 
  ! convert the structure from COO to CSC. Otherwise, we also need
  ! to gather it from the other procs.

  if(transp.eq._qrm_transp) then
     call _qrm_spmat_convert(qrm_mat, graph, 'csr', values=.false., info=err)
     __QRM_INFO_CHECK(err, name,'qrm_spmat_convert',9999)
     graph%irn  => graph%jcn;  nullify(graph%jcn)
     graph%jptr => graph%iptr; nullify(graph%iptr)
     i         = graph%m
     graph%m   = graph%n
     graph%n   = i
     graph%fmt = 'csc'
  else
     call _qrm_spmat_convert(qrm_mat, graph, 'csc', values=.false., info=err)
     __QRM_INFO_CHECK(err, name,'qrm_spmat_convert',9999)
  end if
  call qrm_alloc(dupsmap, graph%m, err)
  __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)
  
  ! make a pass in order to remove duplicates
  dupsmap = 0
  pnt     = 0
  dups    = 0
  savepnt = 1
  do j=1, graph%n
     do ii=graph%jptr(j), graph%jptr(j+1)-1
        i = graph%irn(ii)
        if(dupsmap(i) .eq. j) then
           ! duplicate entry, skip
           dups = dups+1
        else
           ! flag the entry as visited
           dupsmap(i) = j
           pnt = pnt+1
           graph%irn(pnt) = i
        end if
     end do
     graph%jptr(j) = savepnt
     savepnt = pnt+1
  end do
  graph%jptr(graph%n+1)=savepnt
  graph%nz = pnt
  graph%icntl = qrm_mat%icntl
  graph%rcntl = qrm_mat%rcntl

  __QRM_PRNT_DBG('("Number of duplicates in the matrix: ",i10)')dups
  
9999 continue
  ! cleanup and return
  if(err.ne.0) call _qrm_spmat_destroy(graph, all=.true., info=err2)
  if(err2.eq.0) call qrm_dealloc(dupsmap, err2)

  if(present(info)) info = merge(err, err2, err.ne.0)
  return
  
end subroutine _qrm_compute_graph
