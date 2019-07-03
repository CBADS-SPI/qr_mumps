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
!> @file qrm_clean_front.F90
!! This file contains the routines for cleaning of a front
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief This routine performs the cleaning of a front.

!> Cleaning a front means saving the parts corresponding tot he R and
!! Q factors, and then freeing all the memory that is not needed
!! anymore
!!
!! @param[in,out] qrm_mat the whole problem. this obviously contains
!!                        the fornt to be cleaned
!!
!! @param[in,out] front the front to be cleaned
!!
subroutine _qrm_clean_front(qrm_mat, front, info)

  use qrm_mem_mod
  use _qrm_spmat_mod
  use _qrm_fdata_mod
  use _qrm_factorization_mod, protect => _qrm_clean_front  
#if defined (have_starpu)
  use starpu_f_mod
#endif
  implicit none

  type(_qrm_spmat_type), target   :: qrm_mat
  type(_qrm_front_type)           :: front
  integer, optional               :: info
  
  integer                         :: i, j, fc, fr, lc, lr, lt, m, n, tmpmem
  integer(kind=8)                 :: peak, pers
  logical                         :: storeh, storer
  _qrm_data, allocatable          :: tmp(:,:)
  
  ! error management
  integer                         :: err
  character(len=*), parameter     :: name='qrm_clean_front'

  err = 0
  
  if(min(front%m,front%n).le.0) goto 9999

  storer = qrm_mat%icntl(qrm_keeph_).ge.0
  storeh = qrm_mat%icntl(qrm_keeph_).ge.1
  
  ! tmpmem = qrm_tot_mem
  
  ! must use qrm_dealloc without counting
  if(err.eq.0) call qrm_dealloc(front%aiptr, err)
  if(err.eq.0) call qrm_dealloc(front%ajcn,  err)
  if(err.eq.0) call qrm_dealloc(front%aval,  err)
  __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)
  
  if (storer) allocate(front%r((front%npiv-1)/front%mb+1,front%nc))
  if (storeh) allocate(front%h(front%nr,front%nc))
  do j = 1, front%nc
     fc = (j-1)*front%nb+1
     lc = min(j*front%nb, front%n)
     lt = (front%stair(lc)-1)/front%mb+1
     do i = 1, lt
        fr = (i-1)*front%mb+1
        lr = min(fr+size(front%bc(i,j)%c,1)-1, front%m)
        n  = size(front%bc(i,j)%c,2)
        
        ! possible cases:
        ! 1) the block is entirely in R (move_alloc)
        ! 2) the block is entirely in H (move_alloc)
        ! 3) the block is in R, H and/or CB (copies+deallocate)
        ! 4) the block is entirely in CB (deallocate)

        if((lr .lt. min(front%npiv,fc)) .and. storer) then
           ! case 1)
           call qrm_move_alloc(front%bc(i,j)%c, front%r(i,j)%c)
           front%rsize = front%rsize+size(front%bc(i,j)%c)
        else if((fr .gt. lc) .and. storeh) then
           ! case 2)
           call qrm_move_alloc(front%bc(i,j)%c, front%h(i,j)%c)
           front%hsize = front%hsize+size(front%bc(i,j)%c)
        else
           ! case 3)
           ! in this case the lock is necessarily on the diagonal
           if ((fr.le.min(front%npiv,lc)) .and. storer) then
              ! copy the R part
              m = min(lc,front%npiv)-fr+1
              call qrm_alloc(front%r(i,j)%c, m, n, err)
              __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)
              front%r(i,j)%c(1:m,1:n) = front%bc(i,j)%c(1:m,1:n)
              front%rsize = front%rsize+m*n
           end if
           
           if ((lr.ge.fc) .and. storeh) then
              ! copy the H part
              m = lr-fc+1
              call qrm_alloc(front%h(i,j)%c, m, n, err)
              __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)
              front%h(i,j)%c(1:m,1:n) = front%bc(i,j)%c(fc-fr+1:fc-fr+m,1:n)
              front%hsize = front%hsize+m*n
           end if

           call qrm_dealloc(front%bc(i,j)%c)
        end if
        
        if(allocated(front%t(i,j)%c) .and. (size(front%t(i,j)%c,1).ne.front%ib)) then
           ! The size of the T block may be bigger than
           ! necessary. This is due to the trick used to remove false
           ! dependencies between _gemqrt and _tpqrt tasks
           m = front%ib
           call qrm_alloc(tmp, m, size(front%t(i,j)%c,2), err)
           __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)
           tmp(1:m,:) = front%t(i,j)%c(1:m,:)
           call qrm_dealloc(front%t(i,j)%c)
           call qrm_move_alloc(tmp, front%t(i,j)%c)
        end if
        
#if defined(have_starpu)
        if((qrm_mat%adata%small(front%num).eq.0) .and. c_associated(front%bc(i,j)%hdl)) then
           call starpu_f_data_unregister_submit(front%bc(i,j)%hdl)
        end if

        if((qrm_mat%adata%small(front%num).eq.0) .and. c_associated(front%t(i,j)%hdl)) then
           call starpu_f_data_unregister_submit(front%t(i,j)%hdl)
        end if
        if((qrm_mat%adata%small(front%num).eq.0) .and. c_associated(front%t2(i,j)%hdl)) then
           call starpu_f_data_unregister_submit(front%t2(i,j)%hdl)
        end if
#endif
        if(.not.storeh) then
           call qrm_dealloc(front%t(i,j)%c)
           call qrm_dealloc(front%t2(i,j)%c)
        end if

     end do
  end do

  if(.not.storer) then
     call qrm_dealloc(front%cols)
     call qrm_dealloc(front%rows)
     call qrm_dealloc(front%colmap)
     call qrm_dealloc(front%rowmap)
     call qrm_dealloc(front%stair)
  end if

  ! write(*,'(i5,2x,i5," clean -- ",i10,2x,i10,2x,f7.2,2x,i10,2x,i10)')front%num, qrm_mat%adata%small(front%num), &
       ! & tmpmem-qrm_tot_mem, &
       ! & qrm_mat%adata%csize(front%num), &
       ! & 100.d0*real(qrm_mat%adata%csize(front%num)-(tmpmem-qrm_tot_mem))/real(tmpmem-qrm_tot_mem), &
       ! & qrm_tot_mem, qrm_max_mem
  
  if(qrm_mat%adata%small(front%num).ge.0) then
     call qrm_facto_mem_get(qrm_mat%fdata%ma, -qrm_mat%adata%csize(front%num))
  end if

  call qrm_atomic_add(qrm_mat%gstats(qrm_nnz_r_), front%rsize)
  call qrm_atomic_add(qrm_mat%gstats(qrm_nnz_h_), front%hsize)

9999 continue
  if(present(info)) info = err
  return

end subroutine _qrm_clean_front


