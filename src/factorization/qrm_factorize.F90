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
!> @file qrm_factorize.F90
!! This file contains the main factorization driver
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


#include "qrm_common.h"


!> @brief This routine is the main factorization driver
!!
!! @param[in,out] qrm_mat the problem containing the matrix to be factorized.
!!
!! @param[in] transp whether to factorize the input matrix or its
!!            transpose. Accepted values are _qrm_transp or 'n'
!!
subroutine _qrm_factorize(qrm_mat, transp, info)

  use _qrm_spmat_mod
  use qrm_error_mod
  use _qrm_factorization_mod, protect => _qrm_factorize
  use qrm_dscr_mod
  implicit none

  type(_qrm_spmat_type), target   :: qrm_mat
  character, optional, intent(in) :: transp
  integer, optional               :: info

  type(qrm_dscr_type)             :: qrm_dscr

  ! error management
  integer                         :: err
  character(len=*), parameter     :: name='qrm_factorize'

  err = 0

  __QRM_PRNT_DBG('("Entering the factorization driver")')

  ! init the descriptor
  call qrm_dscr_init(qrm_dscr)

  ! initialize the data for the facto
  call _qrm_factorize_async(qrm_dscr, qrm_mat, transp, err)

  ! wait for all the tasks in the dscriptor
  call qrm_barrier(qrm_dscr)
  err = merge(err, qrm_dscr%err_status, err.ne.0)

  call qrm_dscr_destroy(qrm_dscr)

  if(present(info)) info = err
  return

end subroutine _qrm_factorize

!> @brief This routine is the main factorization driver
!!
!! @param[in,out] qrm_mat the problem containing the matrix to be factorized.
!!
!! @param[in] transp whether to factorize the input matrix or its
!!            transpose. Accepted values are _qrm_transp or 'n'
!!
subroutine _qrm_factorize_async(qrm_dscr, qrm_mat, transp, info)

  use _qrm_spmat_mod
  use qrm_error_mod
  use _qrm_factorization_mod, protect => _qrm_factorize_async
  ! use _qrm_fdata_mod
  ! use _qrm_utils_mod
  use qrm_string_mod
  ! use qrm_common_mod
  use qrm_dscr_mod
  implicit none

  type(qrm_dscr_type)             :: qrm_dscr
  type(_qrm_spmat_type), target   :: qrm_mat
  character, optional, intent(in) :: transp
  integer, optional               :: info

  character                       :: itransp
  integer, pointer                :: tmp(:)
  ! error management
  integer                         :: err
  character(len=*), parameter     :: name='qrm_factorize_async'

  err = 0

  ! immediately check if the analysis was done/submitted. Otherwise push an error and return
  if(allocated(qrm_mat%adata)) then
     if(.not. qrm_mat%adata%ok) then
        err = 13
        call qrm_error_print(err, name)
        goto 9999
     end if
  else
     err = 13
     call qrm_error_print(err, name)
     goto 9999
  end if

  ! make sure all the tasks on this matrix (e.g., the analysis) have completed
  call _qrm_sync_spmat(qrm_mat)
  if(qrm_dscr%err_status.ne.0) return
  
  ! check all the control parameters are ok
  call _qrm_check_spmat(qrm_mat, qrm_factorize_, err)
  __QRM_INFO_CHECK(err, name,'qrm_check_spmat',9999)

  if(present(transp)) then
     itransp = qrm_str_tolower(transp)
  else
     itransp = 'n'
  end if

  ! initialize the data for the facto
  call _qrm_factorization_init(qrm_mat, itransp, err)
  __QRM_INFO_CHECK(err, name,'qrm_factorization_init',9999)

  call _qrm_factorization_core(qrm_dscr, qrm_mat, err)
  __QRM_INFO_CHECK(err, name,'qrm_factorization_core',9999)

  ! the factorization was succesfully performed/submitted
  qrm_mat%fdata%ok = .true.

9999 continue
  if(present(info)) info = err
  return

end subroutine _qrm_factorize_async
