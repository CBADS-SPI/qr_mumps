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
!> @file qrm_matnrm.F90
!! this file contains a routine that computes the norm of a sparse matrix
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief This subroutine computes the matrix norm. The return value is
!! a real scalar
!!
!! @param[in] qrm_mat  the inpur A matrix
!! @param[in]   ntype  the norm type. It can be one of these
!!                     1  : 1-norm
!!                     i  : inf-norm
!!                     f  : Frobenius-norm
!! @param[out]    nrm  the output norm
!!
subroutine _qrm_matnrm(qrm_mat, ntype, nrm, info)

  use _qrm_spmat_mod
  use qrm_string_mod
  use qrm_error_mod
  implicit none

  type(_qrm_spmat_type), intent(in) :: qrm_mat
  _qrm_real                         :: nrm
  character                         :: ntype
  integer, optional                         :: info

  _qrm_real, allocatable  :: tmp(:)
  integer                 :: r, c, i
  _qrm_real :: _rxnrm2

  ! error management
  integer                         :: err
  character(len=*), parameter     :: name='qrm_matnrm'

  err = 0

  if(qrm_str_tolower(ntype) .eq. 'i') then
     call qrm_alloc(tmp, qrm_mat%m, err)
     __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)
     tmp = _qrm_zero
     do i=1, qrm_mat%nz
        r = qrm_mat%irn(i)
        tmp(r) = tmp(r)+abs(qrm_mat%val(i))
     end do
     nrm = maxval(tmp)
  else if(qrm_str_tolower(ntype) .eq. '1') then
     call qrm_alloc(tmp, qrm_mat%n, err)
     __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)
     tmp = _qrm_zero
     do i=1, qrm_mat%nz
        c = qrm_mat%jcn(i)
        tmp(c) = tmp(c)+abs(qrm_mat%val(i))
     end do
     nrm = maxval(tmp)
  else if(qrm_str_tolower(ntype) .eq. 'f') then
     nrm = _rxnrm2(qrm_mat%nz, qrm_mat%val, 1)
  else
     err = 15
     call qrm_error_print(err, name)
     goto 9999
  end if

  call qrm_dealloc(tmp)

9999 continue
  if(present(info)) info = err
  return

end subroutine _qrm_matnrm
