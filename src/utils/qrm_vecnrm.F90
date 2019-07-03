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
!> @file qrm_vecnrm.F90
!! This file contains a routine that computes the norm of a vector.
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief This subroutine computes the norm of multiple vectors
!!
!! @param[in]  vec   the input vectors
!!
!! @param[in]  n     the vector size
!!
!! @param[out] nrm   the output norms
!!
!! @param[in]  ntype the norm type. It can be one of these
!!                  1  : 1-norm
!!                  2  : 2-norm
!!                  i  : inf-norm
!!
subroutine _qrm_vecnrm2d(vec, n, ntype, nrm, info)

  use qrm_string_mod
  use qrm_error_mod
  implicit none

  _qrm_data, intent(in)                     :: vec(:,:)
  _qrm_real                                 :: nrm(:)
  integer, intent(in)                       :: n
  character                                 :: ntype
  integer, optional                         :: info

  integer                                   :: i, j
  _qrm_real                                 :: _rxnrm2

  ! error management
  integer                                   :: err
  character(len=*), parameter               :: name='qrm_vecnrm'

  err = 0
  
  nrm = _qrm_rzero

  if(qrm_str_tolower(ntype) .eq. 'i') then
     do j=1, size(vec,2)
        nrm(j) = maxval(abs(vec(:,j)))
     end do
  else if(qrm_str_tolower(ntype) .eq. '1') then
     do j=1, size(vec,2)
        nrm(j) = _qrm_zero
        do i=1, n
           nrm(j) = nrm(j) + abs(vec(i,j))
        end do
     end do
  else if(qrm_str_tolower(ntype) .eq. '2') then
     do j=1, size(vec,2)
        nrm(j) = _rxnrm2(n, vec(1,j), 1)
     end do
  else
     err = 15
     call qrm_error_print(err, name)
  end if
  
9999 continue
  if(present(info)) info = err
  return
  
end subroutine _qrm_vecnrm2d



!> @brief This subroutine computes the norm of a vector
!!
!! @param[in]  vec   the input vector
!!
!! @param[in]  n     the vector size
!!
!! @param[out] nrm   the output norm
!!
!! @param[in]  ntype the norm type. It can be one of these
!!                  1  : 1-norm
!!                  2  : 2-norm
!!                  i  : inf-norm
!!
subroutine _qrm_vecnrm1d(vec, n, ntype, nrm, info)

  use qrm_string_mod
  use qrm_error_mod
  implicit none

  _qrm_data, intent(in)                     :: vec(:)
  _qrm_real                                 :: nrm
  integer, intent(in)                       :: n
  character                                 :: ntype
  integer, optional                         :: info

  integer                                   :: i
  _qrm_real                                 :: _rxnrm2

  ! error management
  integer                                   :: err
  character(len=*), parameter               :: name='qrm_vecnrm'

  err = 0
  
  nrm = _qrm_rzero

  if(qrm_str_tolower(ntype) .eq. 'i') then
     nrm = maxval(abs(vec))
  else if(qrm_str_tolower(ntype) .eq. '1') then
     nrm = _qrm_zero
     do i=1, n
        nrm = nrm + abs(vec(i))
     end do
  else if(qrm_str_tolower(ntype) .eq. '2') then
     nrm = _rxnrm2(n, vec, 1)
  else
     err = 15
     call qrm_error_print(err, name)
  end if
  
9999 continue
  if(present(info)) info = err
  return
  
end subroutine _qrm_vecnrm1d
