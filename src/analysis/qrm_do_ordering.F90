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
!> @file qrm_do_ordering.F90
!! This file contains the routine computes (through different methods) a column permutation
!! of the input matrix in order to reduce fill-in
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


#include "qrm_common.h"


!> @brief This routine computes (through different methods) a column permutation
!! of the input matrix in order to reduce fill-in

!> Supported methods are, currently, COLAMD, SCOTCH and METIS. The user may also
!! provide his own permutation, in which case a check on its validity is done.
!!
!! @param[in] graph    the graph associated to the matrix to be ordered.
!!
!! @param[out] cperm   the new column order
!!
!! @param[in] cperm_in the permutation eventually provided by the user
!!
subroutine _qrm_do_ordering(graph, cperm, cperm_in, info)

  use _qrm_spmat_mod
  use qrm_error_mod
  use _qrm_analysis_mod, savesym => _qrm_do_ordering
  use qrm_mem_mod
  use qrm_const_mod
  use qrm_common_mod
  implicit none

  type(_qrm_spmat_type)       :: graph
  integer                     :: cperm(:)
  integer, pointer            :: cperm_in(:)
  integer, optional           :: info

  integer                     :: iord=0
  integer                     :: i, ts, te, tr

  ! error management
  integer                     :: err
  character(len=*), parameter :: name='qrm_do_ordering'

  err = 0

  call qrm_get(graph, 'qrm_ordering', iord)
  if(iord .eq. qrm_auto_) then
     iord = qrm_choose_ordering()
  end if

  select case(iord)
  case(qrm_natural_)
     ! perm has to be allocated and built explicitly because
     ! we are going to use an equivalent permutation
     do i=1, graph%n
        cperm(i) = i
     end do
  case(qrm_given_)
     ! in this case we just need to check that
     ! the given permutation makes sense
     if(.not. associated(cperm_in)) then
        err = 8; call qrm_error_print(err, name)
        goto 9999
     else
        call qrm_check_cperm(cperm_in, graph%n, err)
        __QRM_INFO_CHECK(err, name,'qrm_check_perm',9999)
        cperm(1:graph%n) = cperm_in(1:graph%n)
     end if
  case(qrm_colamd_)
     ! COLAMD
#if defined (have_colamd)
     call _qrm_do_colamd(graph, cperm, err)
     __QRM_INFO_CHECK(err, name,'qrm_do_colamd',9999)
#else
     err=16; call qrm_error_print(err, name, aed='colamd')
     goto 9999
#endif
  case(qrm_metis_)
     ! METIS
#if defined (have_metis)
     call _qrm_do_metis(graph, cperm, err)
     __QRM_INFO_CHECK(err, name,'qrm_do_metis',9999)
#else
     err=16; call qrm_error_print(err, name, aed='metis')
     goto 9999
#endif
  case(qrm_scotch_)
     ! SCOTCH
#if defined (have_scotch)
     call _qrm_do_scotch(graph, cperm, err)
     __QRM_INFO_CHECK(err, name,'qrm_do_scotch',9999)
#else
     err=16; call qrm_error_print(err, name, aed='scotch')
     goto 9999
#endif
  case default
     err=9; call qrm_error_print(err, name, ied=(/iord/))
     goto 9999
  end select

9999 continue 
  
  if(present(info)) info=err
  return

contains

  function qrm_choose_ordering()
    ! Function: qrm_choose_ordering
    ! This function sets the ordering for the case where an automatic
    ! choice is requested
    !

    integer :: qrm_choose_ordering
    integer :: iord

    iord = 1

#if defined(have_colamd)
    iord = qrm_colamd_
#endif

#if defined(have_scotch)    
    iord = qrm_scotch_
#endif

#if defined(have_metis)    
    iord = qrm_metis_
#endif

    qrm_choose_ordering = iord
    return
    
  end function qrm_choose_ordering


end subroutine _qrm_do_ordering
