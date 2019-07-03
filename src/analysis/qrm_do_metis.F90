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
!> @file qrm_do_metis.F90
!! This file contains the routine that computes a METIS permutation of the input matrix
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> Please refer to:
!!
!! <em>A fast and high quality multilevel scheme for partitioning
!! irregular graphs</em>. George Karypis and Vipin Kumar.  International
!! Conference on Parallel Processing, pp. 113-122, 1995
!! 
!! for the details of the reordering method.
!!
!! @param[in] graph  the graph associated to the matrix to be ordered.
!!
!! @param[out] cperm an integer array containing the new column order
!!
subroutine _qrm_do_metis(graph, cperm, info)

  use _qrm_spmat_mod
  use qrm_error_mod
  use _qrm_analysis_mod, savesym => _qrm_do_metis
  use qrm_mem_mod
  use iso_c_binding

  implicit none

  type(_qrm_spmat_type) :: graph
  integer               :: cperm(:)
  integer, optional     :: info
  
  interface
     subroutine qrm_metis(n, iptr, jcn, cperm, iperm, info) bind(c, name='qrm_metis')
       use iso_c_binding
       integer(c_int)  :: n, info
       integer(c_int)  :: iptr(*), jcn(*), cperm(*), iperm(*)
     end subroutine qrm_metis
  end interface

  integer                     :: i, idx, cnt, tmp, alen, ts, te, tr
  type(_qrm_spmat_type)       :: ata_graph
  integer, allocatable        :: iperm(:)

  ! error management
  integer                     :: err, err2
  character(len=*), parameter :: name='qrm_do_metis'

  err = 0

  call _qrm_ata_graph(graph, ata_graph, err)
  __QRM_INFO_CHECK(err, name,'qrm_ata_graph',9999)

  call qrm_alloc(iperm, graph%n, err)
  __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)

  call qrm_metis(graph%n, ata_graph%iptr, ata_graph%jcn, cperm, iperm, err)
  if(err.ne.0) call qrm_error_print(err, name)
  
9999 continue
  ! cleanup and return
  err2 = 0
  call _qrm_spmat_destroy(ata_graph, all=.true., info=err2)
  if(err2.eq.0) call qrm_dealloc(iperm, err2)

  if(present(info)) info = merge(err, err2, err.ne.0)
  return
  
end subroutine _qrm_do_metis
