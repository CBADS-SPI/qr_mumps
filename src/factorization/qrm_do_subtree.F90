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
!> @file qrm_do_subtree.F90
!! This file contains the routines that take care of factorizing sequentially an entire subtree
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief This subroutine does the sequential factorization of an
!! entire subtree
!!
!! @param[in,out] qrm_mat the whole problem data structure
!!
!! @param[in] root the root of the subtree
!! 
!! @param[in,out] flops a counter for the flops which is updated
!! 
subroutine _qrm_do_subtree(qrm_mat, root, flops, info)
  
  use iso_c_binding
  use qrm_mem_mod
  use _qrm_spmat_mod
  use _qrm_fdata_mod
  use qrm_adata_mod
  use _qrm_utils_mod
  use qrm_common_mod
  use _qrm_factorization_mod, protect => _qrm_do_subtree
  
  implicit none
  
  type(_qrm_spmat_type), target  :: qrm_mat
  type(_qrm_front_type), target  :: root
  real(kind(1.d0))               :: flops
  integer, optional              :: info
  
  type(_qrm_front_type), pointer :: cfront, ffront, front
  type(qrm_adata_type), pointer  :: adata
  type(_qrm_fdata_type), pointer :: fdata
  
  _qrm_data, allocatable         :: rwork(:)
  integer                        :: node, p, c, i, j, cnt, f, m, n, k, ofsa, ib
  integer(kind=8)                :: asize, csize, hsize, rsize
  logical                        :: storer, storeh

  ! error management
  integer                        :: err
  character(len=*), parameter    :: name='qrm_do_subtree'

  err = 0

  ! simplify
  adata => qrm_mat%adata
  fdata => qrm_mat%fdata

  allocate(rwork(adata%rc(root%num)*qrm_mat%icntl(qrm_nb_)))
  
  cnt = 0

  storer = qrm_mat%icntl(qrm_keeph_).ge.0
  storeh = qrm_mat%icntl(qrm_keeph_).ge.1

  ! pointer to the first leaf of the sequential subtree
  node = adata%small(root%num)
  
  subtree: do
     f = adata%torder(node)
     front => fdata%front_list(f)
     ! assemble
     call _qrm_activate_front(qrm_mat, front, info=err)
     __QRM_INFO_CHECK(err, name,'qrm_activate_front',9999)
     call _qrm_init_front(qrm_mat, front, err)
     __QRM_INFO_CHECK(err, name,'qrm_init_front',9999)

     ! factorize
     if(size(rwork).lt.front%n*front%nb) then
        deallocate(rwork)
        allocate(rwork(front%n*front%nb))
     end if

     do j=1, (front%ne-1)/front%nb+1
        ofsa = (j-1)*front%nb+1
        m    = size(front%bc(1,j)%c,1)  - (ofsa-1)
        k    = size(front%bc(1,j)%c,2)
        ib   = min(min(m,k),front%ib)
        call _qrm_geqrt(m, k, ib,                                &
             & front%stair((j-1)*front%nb+1), 0,                 &
             & front%bc(1,j)%c(ofsa,1), size(front%bc(1,j)%c,1), &
             & front%t(1,j)%c(1,1), size(front%t(1,j)%c,1),      &
             & rwork, err)
        k = min(m,k)
        do i=j+1, front%nc
           n = size(front%bc(1,i)%c,2)
           if(n.gt.0) call _qrm_gemqrt('l', _qrm_transp,                    &
                & m, n, k, ib,                                      &
                & front%stair((j-1)*front%nb+1), 0,                 &
                & front%bc(1,j)%c(ofsa,1), size(front%bc(1,j)%c,1), &
                & front%t(1,j)%c(1,1), size(front%t(1,j)%c,1),      &
                & front%bc(1,i)%c(ofsa,1), size(front%bc(1,i)%c,1), &
                & rwork, err)
        end do
     end do

     if(front%num.eq.root%num) exit subtree
     node = node+1

  end do subtree

  ! if root is not the only node in the tree, then we have to release
  ! csize(l) memory, where l is the first leaf of the subtree (see the
  ! comments in qrm_compute_memory.F90
  if(adata%torder(adata%small(root%num)) .ne. root%num) then
     call qrm_facto_mem_get(fdata%ma, -adata%csize(adata%torder(adata%small(root%num))))
  end if

9999 continue ! error management
  if(allocated(rwork)) deallocate(rwork)
  if(present(info)) info = err
  return

end subroutine _qrm_do_subtree







