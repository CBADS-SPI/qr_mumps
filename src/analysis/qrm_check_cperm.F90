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
!> @file qrm_check_cperm.F90
!! This file contains the routine that check whether a provided permutation is good
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief This routine simply checks whether a column permutation provided
!! by the user makes sens.
!! 
!! @param[in] cperm  the permutation to be checked
!! @param[in] n      the size of the permutation verctor, i.e. the number of
!!                   columns in the matrix
!!                   
subroutine qrm_check_cperm(cperm, n, info)
  use qrm_error_mod
  use qrm_mem_mod
  implicit none

  integer                     :: cperm(:)
  integer                     :: n
  integer                     :: info

  integer, allocatable        :: tmp(:)
  integer                     :: i

  ! error management
  integer                     :: err, err2
  character(len=*), parameter :: name='qrm_check_perm'
  
  err = 0; err2 = 0
  
  call qrm_alloc(tmp, n, err)
  __QRM_INFO_CHECK(err, name, 'qrm_alloc', 9999)
  
  tmp = 0

  do i=1, n
     if (cperm(i) .gt. n .or. cperm(i) .lt. 1) then
        err = 8
        goto 9999
     end if
     if(tmp(cperm(i)) .gt. 0) then
        err = 8
        goto 9999
     else
        tmp(cperm(i)) = 1
     end if
  end do

9999 continue
  ! cleanup and return
  call qrm_dealloc(tmp, err2)

  info = merge(err, err2, err.ne.0)
  
  return
  
end subroutine qrm_check_cperm
