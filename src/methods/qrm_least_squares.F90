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
!> @file qrm_analysis_mod.F90
!! This file contains a subroutine for computing the least-squares solution of a problem
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################

#include "qrm_common.h"

!> @brief This routine computes the least-squares solution of a problem

!> This routine computed the least-squares solution of an overdetermined system
!! Ax=b with multiple RHSs.
!!
!! @param[in]     qrm_mat a qrm_spmat_type data which contains the input matrix.
!!                        On output the original data will be unchanged and
!!                        the result of the analysis and factorization phases will
!!                        be stored in the adata and fdata fields, respectively.
!! @param[in,out]       b the RHSs. A 2D array of leading dimension qrm_mat%m. On output
!!                        it will contain Q'*b
!! @param[out]          x the solution, i.e., R\Q'*b
!!
subroutine _qrm_least_squares2d(qrm_mat, b, x, info)
  use _qrm_spmat_mod
  use qrm_error_mod
  use _qrm_analysis_mod
  use _qrm_factorization_mod
  use _qrm_solve_mod
  use _qrm_sdata_mod
  use qrm_dscr_mod
  implicit none

  type(_qrm_spmat_type)            :: qrm_mat
  _qrm_data, intent(inout), target :: b(:,:)
  _qrm_data, intent(out), target   :: x(:,:)
  integer, optional                :: info

  type(_qrm_rhs_type), allocatable :: x_rhs(:), b_rhs(:)
  integer                          :: i, nb, nrhs, n
  type(qrm_dscr_type)              :: qrm_dscr

  ! error management
  integer                          :: err
  character(len=*), parameter      :: name='qrm_least_squares'

  err = 0

  __QRM_PRNT_DBG('("Entering the min-norm driver")')

  call _qrm_check_spmat(qrm_mat, err)
  __QRM_INFO_CHECK(err, name, 'qrm_check_spmat', 9999)

  if(qrm_mat%n .gt. qrm_mat%m) then
     err = 31
     call qrm_error_print(err, name,ied=(/qrm_mat%m,qrm_mat%n/))
     goto 9999
  end if

  ! init the descriptor
  call qrm_dscr_init(qrm_dscr)

  call qrm_analyse_async(qrm_dscr, qrm_mat, 'n', err)
  __QRM_INFO_CHECK(err, name, 'qrm_analyse_async', 9998)
  call qrm_factorize_async(qrm_dscr, qrm_mat, 'n', err)
  __QRM_INFO_CHECK(err, name, 'qrm_factorize_async', 9998)

  ! blocking to deal with multiple rhs
  call qrm_get(qrm_mat, 'qrm_rhsnb', nb)
  nrhs = size(b,2)
  if(nb.le.0) nb = nrhs

  allocate(b_rhs((nrhs-1)/nb+1))
  allocate(x_rhs((nrhs-1)/nb+1))

  do i=1, (nrhs-1)/nb+1
     call _qrm_rhs_init(b_rhs(i), b(:,(i-1)*nb+1:min(i*nb,nrhs)))
     call _qrm_rhs_init(x_rhs(i), x(:,(i-1)*nb+1:min(i*nb,nrhs)))
     call qrm_apply_async(qrm_dscr, qrm_mat, _qrm_transp, b_rhs(i), err)
     __QRM_INFO_CHECK(err, name, 'qrm_apply_async', 9998)
     call qrm_solve_async(qrm_dscr, qrm_mat, 'n', b_rhs(i), x_rhs(i), err)
     __QRM_INFO_CHECK(err, name, 'qrm_solve_async', 9998)
  end do

9998 continue
  
  ! wait for all the tasks in the dscriptor
  call qrm_barrier(qrm_dscr)
  call qrm_dscr_destroy(qrm_dscr)

  err = merge(err, qrm_dscr%err_status, err.ne.0)

  do i=1, (nrhs-1)/nb+1
     call _qrm_rhs_destroy(b_rhs(i))
     call _qrm_rhs_destroy(x_rhs(i))
  end do

9999 continue
  if(present(info)) info = err
  return

end subroutine _qrm_least_squares2d




!> @brief This routine computes the least-squares solution of a problem

!> This routine computed the least-squares solution of an overdetermined system
!! Ax=b with a single RHS.
!!
!! @param[in]     qrm_mat a qrm_spmat_type data which contains the input matrix.
!!                        On output the original data will be unchanged and
!!                        the result of the analysis and factorization phases will
!!                        be stored in the adata and fdata fields, respectively.
!! @param[in,out]       b the RHSs. A 1D array of leading dimension qrm_mat%m. On output
!!                        it will contain Q'*b
!! @param[out]          x the solution, i.e., R\Q'*b
!!
subroutine _qrm_least_squares1d(qrm_mat, b, x, info)
  use _qrm_spmat_mod
  use qrm_error_mod
  use _qrm_utils_mod
  use _qrm_methods_mod, savesym => _qrm_least_squares1d
  implicit none

  type(_qrm_spmat_type)           :: qrm_mat
  _qrm_data, intent(inout)        :: b(:)
  _qrm_data, intent(out)          :: x(:)
  integer, optional               :: info

  _qrm_data, pointer              :: pntb(:,:), pntx(:,:)
  integer                         :: n
  ! error management
  integer                         :: err
  character(len=*), parameter     :: name='qrm_least_squares1d'

  n = size(b,1)
  call _qrm_remap_pnt(b, pntb, n)
  n = size(x,1)
  call _qrm_remap_pnt(x, pntx, n)

  call _qrm_least_squares2d(qrm_mat, pntb, pntx, info)

  return

end subroutine _qrm_least_squares1d
