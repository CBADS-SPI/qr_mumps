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
!> @file qrm_get_r.F90
!! This file contains a subroutine which extracts the R factor
!! resulting from the factorization in sparse form
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################

#include "qrm_common.h"

!> @brief This routine computes the scaled norm of the product A'*r for multiple residuals

!> This routine computes the norm of the scaled product A'*r
!! for a single RHS.
!!
!! @param[in]     qrm_mat a qrm_spmat_type data which contains the R factor to be extracted
!!
!! @param[out]    r  a qrm_spmat_type data which contains the extracted R factor 
!!
subroutine _qrm_get_r(qrm_mat, r, info)
  use qrm_common_mod
  use qrm_error_mod
  use _qrm_fdata_mod
  use _qrm_spmat_mod
  use qrm_mem_mod
  implicit none

  type(_qrm_spmat_type), target  :: qrm_mat
  type(_qrm_spmat_type)          :: r
  integer, optional              :: info

  type(_qrm_front_type), pointer :: front
  integer                        :: cnt, f, i, j, k, m, row, col, rbcnt, rtcnt

  ! error management
  integer                          :: err
  character(len=*), parameter      :: name='qrm_get_r'

  err = 0


  r%nz = qrm_mat%gstats(qrm_nnz_r_)
  r%m = size(qrm_mat%adata%rperm)
  r%n = size(qrm_mat%adata%cperm)

  call qrm_adata_init(r%adata, register=.false.)

  if(err.eq.0) call qrm_alloc(r%irn, r%nz, err)
  if(err.eq.0) call qrm_alloc(r%jcn, r%nz, err)
  if(err.eq.0) call qrm_alloc(r%val, r%nz, err)
  if(err.eq.0) call qrm_alloc(r%adata%rperm, r%m, err)
  if(err.eq.0) call qrm_alloc(r%adata%cperm, r%n, err)
  __QRM_INFO_CHECK(err, name, 'qrm_alloc', 9999)
  r%adata%cperm = qrm_mat%adata%cperm

  rtcnt = 1
  rbcnt = min(r%m,r%n)+1
  cnt   = 1
  do f = 1, qrm_mat%adata%nnodes 
     front => qrm_mat%fdata%front_list(f)
     if(.not.allocated(front%r)) cycle

     r%adata%rperm(rtcnt:rtcnt+front%npiv-1) = front%rows(1:front%npiv)
     r%adata%rperm(rbcnt:rbcnt + front%m-front%ne-1) = front%rows(front%ne+1:front%m)
     rtcnt = rtcnt+front%npiv
     rbcnt = rbcnt + front%m-front%ne
     do col=1, size(front%r,2)
        do row=1,  size(front%r,1)
           if(.not. allocated(front%r(row,col)%c)) cycle
           i = (row-1)*front%mb+1
           j = (col-1)*front%nb+1
           do k=1, size(front%r(row,col)%c, 2)
              m = size(front%r(row,col)%c, 1)
              m = min(m, (j+k-1)-i+1)
              r%irn(cnt:cnt+m-1) = front%rows(i:i+m-1)
              r%jcn(cnt:cnt+m-1) = front%cols(j+k-1)
              r%val(cnt:cnt+m-1) = front%r(row,col)%c(1:m,k)
              cnt=cnt+m
           end do
        end do
     end do
  end do

  r%nz = cnt-1
  if(err.eq.0) call qrm_realloc(r%irn, r%nz, info=err, copy=.true.)
  if(err.eq.0) call qrm_realloc(r%jcn, r%nz, info=err, copy=.true.)
  if(err.eq.0) call qrm_realloc(r%val, r%nz, info=err, copy=.true.)
  __QRM_INFO_CHECK(err, name, 'qrm_realloc', 9999)

  if(rbcnt .ne. r%m+1) then
     __QRM_PRNT_DBG('("_qrm_get_r -- The matrix contains empty rows")')
     r%adata%rperm(rbcnt:r%m) = qrm_mat%adata%rperm(rbcnt:r%m)
  end if

  if(rtcnt.lt.min(r%m,r%n)+1) then
     __QRM_PRNT_ERR('("_qrm_get_r -- The R matrix contains empty rows")')
  end if

  ! if(cnt.ne.r%nz+1) then
  ! __QRM_PRNT_ERR('("_qrm_get_r -- Something went wrong")')
  ! end if

  if(present(info)) info = err
  return
  
9999 continue
  ! cleanup and return
  call qrm_dealloc(r%irn)
  call qrm_dealloc(r%jcn)
  call qrm_dealloc(r%val)
  call qrm_dealloc(r%adata%rperm)
  call qrm_dealloc(r%adata%cperm)
  if(present(info)) info = err

  return
  
end subroutine _qrm_get_r
