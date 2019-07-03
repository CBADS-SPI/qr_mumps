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
!> @file qrm_init_front.F90
!! This file contains a routine which performs the initialization of a front
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


#include "qrm_common.h"


subroutine _qrm_init_front(qrm_mat, front, info)
  use qrm_mem_mod
  use _qrm_spmat_mod
  use _qrm_fdata_mod
  use _qrm_utils_mod
  use qrm_common_mod
  use qrm_sort_mod
  use _qrm_factorization_mod, protect => _qrm_init_front
  implicit none
  
  type(_qrm_spmat_type), target  :: qrm_mat
  type(_qrm_front_type)          :: front
  integer, optional              :: info
  
  integer                        :: i, j, row, col, c
  integer                        :: m, n, npiv, p, ne, bcn, brn, nc, nr, cbcn, mm
  integer                        :: mb, nb, ib

  integer                        :: fnum, ts, ti, taa, tac, tcc, tr
  type(_qrm_front_type), pointer :: cfront
  type(qrm_adata_type), pointer  :: adata
  type(_qrm_fdata_type), pointer :: fdata
  integer(kind=8)                :: peak

  ! error management
  integer                         :: err
  character(len=*), parameter     :: name='qrm_init_front'

  err = 0
  
  cfront => null()

  ! call system_clock(ts)
  
  ! to make things easier
  adata       => qrm_mat%adata
  fdata       => qrm_mat%fdata
  fnum        =  front%num
  front%m     =  adata%nfrows(fnum)
  front%n     =  adata%rc(fnum)
  front%rsize = 0
  front%hsize = 0


  if( (front%n .le. 0) .or. (front%m .le. 0)) then
     ! nothing to do here. Mark the front as done
     front%np = 0
     goto 10
  end if


  ! set up a few variables that come handy in the activation and
  ! factorization of the front
  m          = front%m
  n          = front%n
  nb         = front%nb
  mb         = front%mb
  ib         = front%ib
  ne         = front%ne
  nc         = front%nc
  nr         = front%nr
  npiv       = front%npiv

  cols: do j=1, front%nc
     rows: do i=1, front%nr
        mm = min(front%stair(min(j*nb,n)) - (i-1)*mb, mb)
        if(mm .le. 0) cycle cols

        front%bc(i,j)%c = _qrm_zero
        ! if(allocated(front%t (i,j)%c)) front%t (i,j)%c = _qrm_zero
        ! if(allocated(front%t2(i,j)%c)) front%t2(i,j)%c = _qrm_zero
     end do rows
  end do cols

  ! call system_clock(ti)
  
  
  ! at this point we're ready to assemble the rows from the original matrix
  ! FIXME: deal with the presence of duplicates
  do i=1, front%anrows

     row = front%arowmap(i)
     brn = (row-1)/mb+1
     row = mod(row-1,mb)+1

     !sweep this row and assemble it
     do p=front%aiptr(i), front%aiptr(i+1)-1
        col = front%ajcn(p)
        bcn = (col-1)/nb+1
        col = mod(col-1,nb)+1
        front%bc(brn,bcn)%c(row,col) = front%bc(brn,bcn)%c(row,col)+front%aval(p)
     end do
  end do
  ! call system_clock(taa)

  ! now we can assemble the small children directly
  ! we assume here that because the child is small it only has one block
  do p = adata%childptr(fnum), adata%childptr(fnum+1)-1
     c = adata%child(p)
     cfront => fdata%front_list(c)

     ! ne is the number of Householder vectors computed on the
     ! child c. npiv is the number of fully assembled pivots in c
     ne   = cfront%ne
     npiv = cfront%npiv
     if(npiv.eq.ne) cycle

     ! fill-up front%rows using the rowmap computed in the activation
     do i=npiv+1, ne
        row = cfront%rowmap(i-npiv)
        front%rows(row) = cfront%rows(i)
     end do
     
     if(adata%small(c) .ne. 0) then

        do j=npiv+1, cfront%n
           ! this is the column mapping on the child. cfront%colmap(k)=j
           ! means that the k-th column of cfront will be assembled into the
           ! j-th column of front

           col    = cfront%colmap(j-npiv)
           bcn    = (col-1)/nb+1
           col    = mod(col-1,nb)+1
           cbcn   = (j-1)/cfront%nb+1
           
           ! fill in the front with coefficients front child front cfront
           do i=npiv+1, min(j,ne)
              row = cfront%rowmap(i-npiv)
              brn = (row-1)/mb+1
              row = mod(row-1,mb)+1
              front%bc(brn,bcn)%c(row,col) = cfront%bc(1,cbcn)%c(i,j-(cbcn-1)*cfront%nb)
           end do
        end do
     end if
  end do

  ! call system_clock(tac)

  
10 continue


  do p = adata%childptr(fnum), adata%childptr(fnum+1)-1
     c = adata%child(p)
     if(adata%small(c) .lt. 0) then
        cfront => fdata%front_list(c)
        peak = qrm_tot_mem
        call _qrm_clean_front(qrm_mat, cfront, err)
        __QRM_INFO_CHECK(err, name, 'qrm_clean_front', 9999)
     end if
  end do

  ! call system_clock(tcc,tr)


  ! if(adata%small(fnum) .eq. 0) then
  !    write(*,'(i6,"  Init times -- total:",f8.5,"  zero:",f8.5,"  aasm:",f8.5,"  casm:",f8.5,"  ccln:",f8.5)') &
  !         & front%num, &
  !         & real(tcc-ts)/real(tr), &
  !         & real(ti-ts)/real(tr), &
  !         & real(taa-ti)/real(tr), &
  !         & real(tac-taa)/real(tr), &
  !         & real(tcc-tac)/real(tr)
  ! end if
  
9999 continue
  if(present(info)) info = err
  return

end subroutine _qrm_init_front
