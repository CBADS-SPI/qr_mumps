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

#include "qrm_common.h"

subroutine _qrm_compute_memory(qrm_mat, transp, info)

  use _qrm_spmat_mod
  use qrm_adata_mod
  use qrm_error_mod
  use qrm_mem_mod
  use qrm_memhandling_mod
  use qrm_common_mod
  use _qrm_utils_mod
  implicit none
  
  type(_qrm_spmat_type), target :: qrm_mat
  character                     :: transp
  integer, optional             :: info
  
  type(qrm_adata_type), pointer :: adata
  integer(kind=8)               :: peak, cons, subtree_peak, subtree_res, subtree_cons, cfree
  integer                       :: i, f, c, j, pfront
  integer(kind=8), allocatable  :: tmp(:)
  integer                       :: mb, nb, ib, bh, ne, np
  integer(kind=8)               :: frsize, fhsize, size_arrowheads
  logical                       :: storeh, storer

  ! error management
  integer                       :: err
  character(len=*), parameter   :: name='qrm_compute_memory'

  err = 0
  
  ! just to simplify
  adata    => qrm_mat%adata

  ! determine mb, nb, ib and bh
  call qrm_get(qrm_mat, 'qrm_mb', mb)
  call qrm_get(qrm_mat, 'qrm_nb', nb)
  call qrm_get(qrm_mat, 'qrm_ib', ib)
  call qrm_get(qrm_mat, 'qrm_bh', bh)
  storer = qrm_mat%icntl(qrm_keeph_).ge.0
  storeh = qrm_mat%icntl(qrm_keeph_).ge.1
  
  ! update the memory footprint of small nodes
  ! do f=1, adata%nnodes
     ! if(adata%small(f).ne.0) then
        ! ne = min(adata%rc(f),adata%nfrows(f))
        ! np = adata%cp_ptr(f+1)-adata%cp_ptr(f)
        ! call qrm_get_front_mem(adata%nfrows(f), adata%rc(f), &
             ! & adata%nfrows(f), adata%rc(f), ib, bh, np, &
             ! & adata%asize(f), adata%csize(f), fhsize, frsize, &
             ! & storer, storeh)
     ! end if
  ! end do
  
  cons         = 0
  ! peak is initialized to the current consumption. this is ok because
  ! this routine is the last to be called in te analsys phase
  peak         = 0
  pfront       = 0
  subtree_peak = 0
  subtree_res  = 0
  subtree_cons = 0

  ! add the arrowheads. Note that the arrowheads are freed in qrm_clean
  ! front but not counted in cmem which leads to an overestimation of
  ! the mem consumption (there are also other things that contribute
  ! such as the T matrices)

  if(transp.eq.'n') then
     size_arrowheads = 2*qrm_mat%m*_qrm_sizeof_i ! aiptr + arowmap
  else if (transp.eq._qrm_transp) then
     size_arrowheads = 2*qrm_mat%n*_qrm_sizeof_i ! aiptr + arowmap
  end if
  size_arrowheads = size_arrowheads +   &
       & (adata%nnodes +                &        ! one extra for aiptr
       & qrm_mat%nz)                    &        ! ajcn
       & *_qrm_sizeof_i +               &    
       & qrm_mat%nz*_qrm_sizeof_data             ! aval

  cons = cons + size_arrowheads

  do i=1, adata%nnodes
     f = adata%torder(i)
     cons = cons+adata%asize(f)
     if(cons.gt.peak) then
        peak = cons
        pfront = f
     end if
     ! write(*,'(i5," activ -- ",i10,2x,i10,2x,i10,2x,i10)')f, cons, peak, adata%asize(f), adata%csize(f)

     cfree=0
     do j=adata%childptr(f), adata%childptr(f+1)-1
        c = qrm_mat%adata%child(j)
        cfree = cfree+adata%csize(c)
        cons = cons - adata%csize(c)
        ! write(*,'(i5," clean -- ",i10,2x,i10,2x,i10)')c, adata%csize(c),cons, peak
     end do

     ! Here's the logic of asize and csize for subtrees:
     !
     ! if f is the root of a sequential subtree but is not the only
     ! node in a subtree, asize(f) is the peak consumption of the
     ! subtree and csize(l), where l is the first leaf of the subtree,
     ! is the amount of memory freed (wrt the peak) when all the nodes
     ! in the tree are deactivated except the root node. Therefore, in
     ! order to process the subtree we have to reserve asize(f); when
     ! we get to the root of the subtree we free csize(l). When the
     ! root is assembled into its father we free csize(f).
     

     if(adata%small(f).ne.0) then
        ! we are inside a sequential subtree. update it's memory
        ! footprint
        subtree_cons                                    = subtree_cons + adata%asize(f)
        if(subtree_cons .gt. subtree_peak) subtree_peak = subtree_cons
        subtree_cons                                    = subtree_cons-cfree
        subtree_res                                     = subtree_res + adata%asize(f)-adata%csize(f)
        if(adata%small(f).gt.0) then
           ! we are on the root of a sequential subtree.
           if (adata%small(f).ne.i) then
              ! replace its asize and csize with those of the whole
              ! subtree
              adata%asize(f) = subtree_peak
              ! csize(l) is equal to the peak - whathever is left when
              ! we have assembled the root and feed its children
              adata%csize(adata%torder(adata%small(f))) = subtree_peak - subtree_res - adata%csize(f)
           end if
           subtree_peak   = 0
           subtree_res    = 0
           subtree_cons   = 0
        end if
     end if
     
  end do

  qrm_mat%gstats(qrm_e_facto_mempeak_) = peak
  ! write(*,'("Estimated memory peak (MB): ",i12,"   reached on front:",i6)')peak, pfront

9999 continue

  if(present(info)) info = err
  return

end subroutine _qrm_compute_memory
