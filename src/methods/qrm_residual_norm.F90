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
!> @file qrm_residual_norm.F90
!! This file contains a subroutine for computing the scaled norm of the residual
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################

#include "qrm_common.h"

!> @brief This routine computes the scaled norm of multiple residuals

!> This routine computes the norm of the scaled residual, i.e., ||r||/||A||
!! for a single RHS.
!!
!! @param[in]     qrm_mat a qrm_spmat_type data which contains the input matrix.
!!
!! @param[in,out]       b the RHSs. A 1D array of leading dimension qrm_mat%m. On output
!!                        it will contain the residual
!!
!! @param[in]           x the solution
!!
!! @param[out]        nrm the residual norm
!!
subroutine _qrm_residual_norm2d(qrm_mat, b, x, nrm, info)
  use _qrm_spmat_mod
  use qrm_error_mod
  use _qrm_utils_mod
  implicit none

  type(_qrm_spmat_type)  :: qrm_mat
  _qrm_data              :: b(:,:), x(:,:)
  _qrm_real              :: nrm(:)
  integer, optional      :: info

  _qrm_real, allocatable :: nrmb(:), nrmx(:)
  _qrm_real              :: nrma
  integer                :: nrhs

  ! error management
  integer                         :: err
  character(len=*), parameter     :: name='qrm_residual_norm'

  err = 0

  call _qrm_check_spmat(qrm_mat, err)
  __QRM_INFO_CHECK(err, name,'qrm_check_spmat',9999)

  nrhs = min(size(b,2),size(x,2))

  call qrm_alloc(nrmb, nrhs, err)
  call qrm_alloc(nrmx, nrhs, err)
  __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)

  call qrm_vecnrm(b(:,1:nrhs), qrm_mat%m, 'i', nrmb)
  call qrm_vecnrm(x(:,1:nrhs), qrm_mat%n, 'i', nrmx)

  ! compute the residual
  call qrm_matmul(qrm_mat, 'n', -_qrm_one, x, _qrm_one, b)
  call qrm_matnrm(qrm_mat, 'i', nrma)
  call qrm_vecnrm(b, qrm_mat%m, 'i', nrm)

  nrmb = nrmb + nrma*nrmx

  nrm = nrm/nrmb

  call qrm_dealloc(nrmx)
  call qrm_dealloc(nrmb)

9999 continue ! error management
  if(present(info)) info = err
  return

end subroutine _qrm_residual_norm2d



!> @brief This routine computes the scaled norm of the residual

!> This routine computes the norm of the scaled residual, i.e., ||r||/||A||
!! for a single RHS.
!!
!! @param[in]     qrm_mat a qrm_spmat_type data which contains the input matrix.
!!
!! @param[in,out]       b the RHSs. A 1D array of leading dimension qrm_mat%m. On output
!!                        it will contain the residual
!!
!! @param[in]           x the solution
!!
!! @param[out]        nrm the residual norm
!!
subroutine _qrm_residual_norm1d(qrm_mat, b, x, nrm, info)
  use _qrm_spmat_mod
  use qrm_error_mod
  use _qrm_utils_mod
  implicit none

  type(_qrm_spmat_type) :: qrm_mat
  _qrm_data             :: b(:), x(:)
  _qrm_real             :: nrm
  integer, optional     :: info

  _qrm_real             :: nrmb, nrmx
  _qrm_real             :: nrma


  ! error management
  integer                         :: err
  character(len=*), parameter     :: name='qrm_residual_norm'

  err = 0

  call _qrm_check_spmat(qrm_mat, err)
  __QRM_INFO_CHECK(err, name, 'qrm_check_spmat', 9999)

  call qrm_vecnrm(b, qrm_mat%m, 'i', nrmb)
  call qrm_vecnrm(x, qrm_mat%n, 'i', nrmx)

  ! compute the residual
  call qrm_matmul(qrm_mat, 'n', -_qrm_one, x, _qrm_one, b)
  call qrm_matnrm(qrm_mat, 'i', nrma)
  call qrm_vecnrm(b, qrm_mat%m, 'i', nrm)

  nrmb = nrmb + nrma*nrmx
  nrm = nrm/nrmb

9999 continue ! error management
  if(present(info)) info = err
  return

end subroutine _qrm_residual_norm1d
