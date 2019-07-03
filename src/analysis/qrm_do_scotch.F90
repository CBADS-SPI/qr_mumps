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
!> @file qrm_do_scotch.F90
!! This file contains the routine that computes a SCOTCH permutation of the input matrix
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
!! <em>"SCOTCH 5.1 User's guide</em>.
!! Technical report, LaBRI, September 2008.
!! F. Pellegrini.
!!
!! for the details of the reordering method.
!!
!! @param[in] graph  the graph associated to the matrix to be ordered.
!!
!! @param[out] cperm an integer array containing the new column order
!!
subroutine _qrm_do_scotch(graph, cperm, info)

  use _qrm_spmat_mod
  use qrm_error_mod
  use _qrm_analysis_mod, savesym => _qrm_do_scotch
  use qrm_mem_mod
  use iso_c_binding
  implicit none

#if defined(have_scotch)
  include "scotchf.h"
#endif

  type(_qrm_spmat_type)       :: graph
  integer                     :: cperm(:)
  integer, optional           :: info

  ! TODO: fix this with conditional compilation
#if defined(have_scotch)
  integer                     :: i, cblknbr, ts, te, tr
  type(_qrm_spmat_type)       :: ata_graph
  integer, allocatable        :: iperm(:)
  real(kind(1.d0))            :: grafdat(scotch_graphdim), stradat(scotch_stratdim), orderdat(scotch_orderdim)

  ! error management
  integer                     :: err
  character(len=*), parameter :: name='qrm_do_scotch'

  err = 0

  call _qrm_ata_graph(graph, ata_graph, err)
  __QRM_INFO_CHECK(err, name,'qrm_ata_graph',9999)
  
  info = 0
  call scotchfgraphinit(grafdat, err)
  call scotchfstratinit(stradat, err)
  if(err .ne. 0) then
     err = 19
     call qrm_error_print(err, name)
     goto 9999
  end if

  call scotchfgraphbuild(grafdat, 1, ata_graph%n, ata_graph%iptr(1), &
       & ata_graph%iptr(2), ata_graph%iptr, ata_graph%iptr, ata_graph%nz, &
       & ata_graph%jcn, ata_graph%jcn, err)
  if(err .ne. 0) then
     err = 19
     call qrm_error_print(err, name)
     goto 9999
  end if

  call scotchfgraphorder(grafdat, stradat, grafdat, cperm, cblknbr, &
       & grafdat, grafdat, info)
  if(err .ne. 0) then
     err = 19
     call qrm_error_print(err, name)
     goto 9999
  end if

  call scotchfgraphexit(grafdat)
  call scotchfstratexit(stradat)

9999 continue ! error management

  call _qrm_spmat_destroy(ata_graph, all=.true.)

  if(present(info)) info=err
  return

#endif

  return
end subroutine _qrm_do_scotch
