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
!> @file qrm_solve_tasks.F90
!! This file contains all the kernel for the solve phase
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################

#include "qrm_common.h"

!==========================================================================================
!==========================================================================================
subroutine _qrm_front_qt(front, b, work)
  use _qrm_fdata_mod
  use _qrm_sdata_mod
  implicit none

  type(_qrm_front_type)           :: front
  _qrm_data                       :: work(:,:)
  type(_qrm_rhs_type)             :: b

  integer                         :: lt, ft, lct, nsteps, s, bi, bk, bl
  integer                         :: i, m, n, k, ib, l, node, info

  ! error management
  character(len=*), parameter     :: name='front_qt'

  ! shortcut
  if (min(front%m, front%n) .le. 0) goto 9999

  node = front%num
  
  n = size(b%p,2)

  ! write(*,'(" =--> Apply Q''  : ",i4)')front%num

  lt = (min(front%m,front%n)-1)/front%nb+1
  do bk=1, lt
     ft = ((bk-1)*front%nb)/front%mb+1
     lct = (front%stair(min(bk*front%nb, front%n))-1)/front%mb+1

     do bi = ft, lct, front%bh
        i  = max((bk-1)*front%nb+1 - (bi-1)*front%mb,1)
        m  = size(front%h(bi,bk)%c,1)
        k  = size(front%h(bi,bk)%c,2)
        ib = min(front%ib,k)
        ib = min(ib,m)
        call _qrm_gemqrt('l', _qrm_transp,                                  &
             & m, n, k, ib,                                         &
             & front%stair((bk-1)*front%nb+1), (bi-1)*front%mb+i-1, &
             & front%h(bi,bk)%c(1,1), size(front%h(bi,bk)%c,1),     &
             & front%t(bi,bk)%c(1,1), size(front%t(bi,bk)%c,1),     &
             & b%front_rhs(node)%c((bi-1)*front%mb+i, 1), front%m,  &
             & work, info)
        do bl=bi+1, min(bi+front%bh-1,lct)
           m  = size(front%h(bl,bk)%c,1)
           k  = size(front%h(bi,bk)%c,2)
           l  = 0
           ib = min(front%ib,k)

           call _qrm_tpmqrt('l', _qrm_transp,                                  &
                & m, n, k, l, ib,                                      &
                & front%stair((bk-1)*front%nb+1), (bi-1)*front%mb+i-1, &
                & (bl-1)*front%mb,                                     &
                & front%h(bl,bk)%c(1,1), size(front%h(bl,bk)%c,1),     &
                & front%t2(bl,bk)%c(1,1), size(front%t2(bl,bk)%c,1),   &
                & b%front_rhs(node)%c((bi-1)*front%mb+i, 1), front%m,  &
                & b%front_rhs(node)%c((bl-1)*front%mb+1, 1), front%m,  &
                & work, info)
        end do
     end do
     s = front%bh
     do while (s.le.lct-ft+1)
        do bi = ft, lct-s, 2*s
           i  = max((bk-1)*front%nb+1 - (bi-1)*front%mb,1)
           bl = bi+s
           if(bl.le.lct) then
              m  = min(size(front%h(bl,bk)%c,1), size(front%h(bl,bk)%c,2))
              k  = size(front%h(bi,bk)%c,2)
              l  = m
              ib = min(front%ib,k)

              call _qrm_tpmqrt('l', _qrm_transp, &
                   & m, n, k, l, ib, &
                   & front%stair((bk-1)*front%nb+1), (bi-1)*front%mb+i-1, (bl-1)*front%mb, &
                   & front%h(bl,bk)%c(1,1), size(front%h(bl,bk)%c,1), &
                   & front%t2(bl,bk)%c(1,1), size(front%t2(bl,bk)%c,1), &
                   & b%front_rhs(node)%c((bi-1)*front%mb+i, 1), front%m, &
                   & b%front_rhs(node)%c((bl-1)*front%mb+1, 1), front%m, &
                   & work, info)
           end if
        end do
        s = s*2
     end do
  end do

  ! scatter front%b into b
  b%p(front%rows(1:front%npiv),:) = b%front_rhs(node)%c(1:front%npiv,:)
  if(front%m.gt.front%ne) then
     b%p(front%rows(front%ne+1:front%m),:) = b%front_rhs(node)%c(front%ne+1:front%m,:)
  end if

9999 continue

  return 

end subroutine _qrm_front_qt


subroutine _qrm_assemble_qt(qrm_mat, front, b, info)
  use _qrm_spmat_mod
  use qrm_mem_mod
  use qrm_error_mod
  use _qrm_sdata_mod
  use _qrm_fdata_mod
  implicit none

  type(_qrm_spmat_type), target  :: qrm_mat
  type(_qrm_front_type)          :: front
  type(_qrm_rhs_type)            :: b
  integer, optional              :: info

  type(qrm_adata_type), pointer  :: adata
  type(_qrm_front_type), pointer :: child
  integer                        :: n, c, node

  ! error management
  integer                        :: err
  character(len=*), parameter    :: name='qrm_assemble_qt'

  err = 0

  node = front%num
  
  adata => qrm_mat%adata
  if(front%ne .gt. 0) then
     n = size(b%p,2)
     call qrm_alloc(b%front_rhs(node)%c, front%m, n, err)
     __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)
     if(front%anrows.gt.0)b%front_rhs(node)%c(front%arowmap,:) = b%p(front%rows(front%arowmap),:)
  end if

  do c = adata%childptr(node), adata%childptr(node+1)-1
     child => qrm_mat%fdata%front_list(adata%child(c))

     ! the double condition below should be redundant but you never know...
     if ((child%ne.gt.child%npiv) .and. (front%ne.gt.0)) &
          & b%front_rhs(node)%c(child%rowmap,:) = &
          & b%front_rhs(child%num)%c(child%npiv+1:child%ne,:)
     call qrm_dealloc(b%front_rhs(child%num)%c)
  end do


9999 continue

  if(present(info)) info = err
  return

end subroutine _qrm_assemble_qt





!==========================================================================================
!==========================================================================================
subroutine _qrm_front_r(front, b, x)

  use _qrm_fdata_mod
  use _qrm_sdata_mod
  implicit none

  type(_qrm_front_type)           :: front
  type(_qrm_rhs_type)             :: b, x

  integer                         :: bi, bj, bk, npi, npj
  integer                         :: m, n, k, i, j, node, info

  ! shortcut
  if (min(front%m, front%n) .le.0) goto 9999
  if (front%npiv .le. 0) goto 9999

  node = front%num
  
  ! write(*,'(" =--> Solve R  : ",i4)') front%num

  ! 1st part, multiply by the blocks in npiv+1:n
  n   = size(x%p,2)
  npj = (front%npiv)/front%nb + 1
  npi = (front%npiv-1)/front%mb + 1
  do bj=npj, front%nc
     do bi=1, npi
        m = min(bi*front%mb,front%npiv) - (bi-1)*front%mb  ! size(front%r(bi,bj)%c,1)
        j = max(1,front%npiv-(bj-1)*front%nb+1)

        k = size(front%r(bi,bj)%c,2)-j+1
        if(k.gt.0) call _xgemm('n', 'n', &
             & m, n, k, &
             & _qrm_mone, &
             & front%r(bi,bj)%c(1,j), size(front%r(bi,bj)%c,1), &
             & x%front_rhs(node)%c((bj-1)*front%nb+j,1), front%n, &
             & _qrm_one, &
             & x%front_rhs(node)%c((bi-1)*front%mb+1,1), front%n)
     end do
  end do


  ! 2nd part, do the trsm with blocks in 1:npiv
  npj = (front%npiv-1)/front%nb + 1
  do bj=npj, 1, -1
     i  = (bj-1)*front%nb+1
     bi = (i-1)/front%mb +1
     m  = min(front%npiv-i+1,size(front%r(bi,bj)%c,2))
     i  = i - (bi-1)*front%mb
     call _xtrsm('l', 'u', 'n', 'n', &
          & m, n, &
          & _qrm_one, &
          & front%r(bi,bj)%c(i,1), size(front%r(bi,bj)%c,1), &
          & x%front_rhs(node)%c((bj-1)*front%nb+1,1), front%n)

     k = m
     do while (bi.gt.0)
        m = min((bj-1)*front%nb - (bi-1)*front%mb, size(front%r(bi,bj)%c,1))
        if(m.gt.0) call _xgemm('n', 'n', &
             & m, n, k, &
             & _qrm_mone, &
             & front%r(bi,bj)%c(1,1), size(front%r(bi,bj)%c,1), &
             & x%front_rhs(node)%c((bj-1)*front%nb+1,1), front%n, &
             & _qrm_one, &
             & x%front_rhs(node)%c((bi-1)*front%mb+1,1), front%n)
        bi = bi-1
     end do
  end do


9999 continue 
  return

end subroutine _qrm_front_r


subroutine _qrm_assemble_r(qrm_mat, front, b, x, info)
  use _qrm_spmat_mod
  use qrm_mem_mod
  use qrm_error_mod
  use _qrm_fdata_mod
  use _qrm_sdata_mod
  implicit none

  type(_qrm_spmat_type), target  :: qrm_mat
  type(_qrm_front_type)          :: front
  type(_qrm_rhs_type)            :: b, x
  integer, optional              :: info

  type(qrm_adata_type), pointer  :: adata
  type(_qrm_front_type), pointer :: child
  integer                        :: n, c, node

  ! error management
  integer                        :: err
  character(len=*), parameter    :: name='qrm_assemble_r'

  err = 0

  adata => qrm_mat%adata

  node = front%num
  n = size(b%p,2)

  do c = adata%childptr(node), adata%childptr(node+1)-1
     child => qrm_mat%fdata%front_list(adata%child(c))
     call qrm_alloc(x%front_rhs(child%num)%c, child%n, n, err)
     __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)

     x%front_rhs(child%num)%c(1:child%npiv,:) = b%p(child%rows(1:child%npiv),:)
     if(child%npiv.lt.child%n) x%front_rhs(child%num)%c(child%npiv+1:child%n,:) = x%front_rhs(node)%c(child%colmap,:)

  end do

  if(front%ne .le. 0) goto 9999

  ! scatter front%rhs into x
  x%p(front%cols(1:front%npiv),:) = x%front_rhs(node)%c(1:front%npiv,:)

  call qrm_dealloc(x%front_rhs(node)%c)

9999 continue

  if(present(info)) info = err
  return

end subroutine _qrm_assemble_r



!==========================================================================================
!==========================================================================================
subroutine _qrm_front_q(front, b, work)
  use _qrm_utils_mod
  use _qrm_fdata_mod
  use qrm_mem_mod
  use qrm_common_mod
  use _qrm_sdata_mod
  implicit none

  type(_qrm_front_type)           :: front
  type(_qrm_rhs_type)             :: b
  _qrm_data                       :: work(:,:)

  integer                         :: lt, ft, lct, nsteps, s, bi, bk, bl
  integer                         :: i, m, n, k, f, ib, l, node, info

  ! write(*,'(" =--> Apply Q  : ",i4)')front%num

  ! shortcut
  if (min(front%m, front%n) .le. 0) goto 9999

  node = front%num
  n = size(b%p,2)

  lt = (min(front%m,front%n)-1)/front%nb+1
  do bk=lt, 1, -1
     ft = ((bk-1)*front%nb)/front%mb+1
     lct = (front%stair(min(bk*front%nb, front%n))-1)/front%mb+1
     s = front%bh
     do while (s.le.lct-ft+1)
        s = s*2
     end do
     do while (s.ge.front%bh)
        do bi = ft, lct-s, 2*s
           i  = max((bk-1)*front%nb+1 - (bi-1)*front%mb,1)
           bl = bi+s
           if(bl.le.lct) then
              m  = min(size(front%h(bl,bk)%c,1), size(front%h(bl,bk)%c,2))
              k  = size(front%h(bi,bk)%c,2)
              l  = m
              ib = min(front%ib,k)

              call _qrm_tpmqrt('l', 'n', &
                   & m, n, k, l, ib, &
                   & front%stair((bk-1)*front%nb+1), (bi-1)*front%mb+i-1, (bl-1)*front%mb, &
                   & front%h(bl,bk)%c(1,1), size(front%h(bl,bk)%c,1), &
                   & front%t2(bl,bk)%c(1,1), size(front%t2(bl,bk)%c,1), &
                   & b%front_rhs(node)%c((bi-1)*front%mb+i, 1), front%m, &
                   & b%front_rhs(node)%c((bl-1)*front%mb+1, 1), front%m, &
                   & work, info)
           end if
        end do
        s = s/2
     end do

     do bi = ft, lct, front%bh
        i  = max((bk-1)*front%nb+1 - (bi-1)*front%mb,1)

        do bl=min(bi+front%bh-1,lct), bi+1, -1
           m  = size(front%h(bl,bk)%c,1)
           k  = size(front%h(bi,bk)%c,2)
           l  = 0
           ib = min(front%ib,k)
           call _qrm_tpmqrt('l', 'n', &
                & m, n, k, l, ib, &
                & front%stair((bk-1)*front%nb+1), (bi-1)*front%mb+i-1, (bl-1)*front%mb, &
                & front%h(bl,bk)%c(1,1), size(front%h(bl,bk)%c,1), &
                & front%t2(bl,bk)%c(1,1), size(front%t2(bl,bk)%c,1), &
                & b%front_rhs(node)%c((bi-1)*front%mb+i, 1), front%m, &
                & b%front_rhs(node)%c((bl-1)*front%mb+1, 1), front%m, &
                & work, info)
        end do


        m  = size(front%h(bi,bk)%c,1)
        k  = size(front%h(bi,bk)%c,2)
        ib = min(front%ib,k)
        ib = min(ib,m)
        call _qrm_gemqrt('l', 'n', &
             & m, n, k, ib, &
             & front%stair((bk-1)*front%nb+1), (bi-1)*front%mb+i-1, &
             & front%h(bi,bk)%c(1,1), size(front%h(bi,bk)%c,1), &
             & front%t(bi,bk)%c(1,1), size(front%t(bi,bk)%c,1), &
             & b%front_rhs(node)%c((bi-1)*front%mb+i, 1), front%m, &
             & work, info)
     end do

  end do

9999 continue
  return 

end subroutine _qrm_front_q

subroutine _qrm_assemble_q(qrm_mat, front, b, info)
  use _qrm_spmat_mod
  use qrm_mem_mod
  use qrm_error_mod
  use qrm_common_mod
  use _qrm_sdata_mod
  implicit none

  type(_qrm_spmat_type), target  :: qrm_mat
  type(_qrm_front_type)          :: front
  type(_qrm_rhs_type)            :: b
  integer, optional              :: info

  type(qrm_adata_type), pointer  :: adata
  type(_qrm_front_type), pointer :: child
  integer                        :: n, c, node

  ! error management
  integer                        :: err
  character(len=*), parameter    :: name='qrm_assemble_q'

  err = 0

  adata => qrm_mat%adata

  n = size(b%p,2)
  node = front%num
  
  do c = adata%childptr(front%num), adata%childptr(front%num+1)-1
     child => qrm_mat%fdata%front_list(adata%child(c))
     call qrm_alloc(b%front_rhs(child%num)%c, child%m, n, err)
     __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)

     b%front_rhs(child%num)%c(1:child%npiv,:) = b%p(child%rows(1:child%npiv),:)
     if(child%m.ge.child%ne) b%front_rhs(child%num)%c(child%ne+1:child%m,:) = b%p(child%rows(child%ne+1:child%m),:)
     if(child%ne.gt.child%npiv) b%front_rhs(child%num)%c(child%npiv+1:child%ne,:) = b%front_rhs(node)%c(child%rowmap,:)

  end do

  if(front%ne .le. 0) goto 9999

  ! scatter front%b into b
  b%p(front%rows,:) = b%front_rhs(node)%c

  call qrm_dealloc(b%front_rhs(node)%c)

9999 continue
  
  if(present(info)) info = err
  return

end subroutine _qrm_assemble_q



!==========================================================================================
!==========================================================================================
subroutine _qrm_front_rt(front, b, x)
  use _qrm_fdata_mod
  use _qrm_sdata_mod
  implicit none

  type(_qrm_front_type)           :: front
  type(_qrm_rhs_type)             :: b, x

#if defined (sprec) || defined (dprec)
  character, parameter            :: tr=_qrm_transp, notr='n'
#elif defined (cprec) || defined (zprec)
  character, parameter            :: tr='c', notr='n'
#endif

  integer                         :: bi, bj, bk, npi, npj, m, n, k, i, j, node, info

  ! shortcut
  if (min(front%m, front%n) .le. 0) goto 10
  if (front%npiv .le. 0) goto 10

  node = front%num
  ! write(*,'(" =--> Solve RT  : ",i4)') front%num

  n = size(x%p,2)

  ! 1st part, do the trsm with blocks in 1:npiv
  npj = (front%npiv-1)/front%nb + 1
  do bj=1, npj
     i  = (bj-1)*front%nb+1
     m  = min(front%npiv-i+1,size(front%r(1,bj)%c,2))

     do bi=1, (i-1)/front%mb+1
        k = min((bj-1)*front%nb - (bi-1)*front%mb, size(front%r(bi,bj)%c,1))
        if(k.gt.0) call _xgemm(tr, 'n', &
             & m, n, k, &
             & _qrm_mone, &
             & front%r(bi,bj)%c(1,1), size(front%r(bi,bj)%c,1), &
             & x%front_rhs(node)%c((bi-1)*front%mb+1,1), front%n, &
             & _qrm_one, &
             & x%front_rhs(node)%c((bj-1)*front%nb+1,1), front%n)
     end do
     bi = (i-1)/front%mb +1
     i  = i - (bi-1)*front%mb
     call _xtrsm('l', 'u', tr, 'n', &
          & m, n, &
          & _qrm_one, &
          & front%r(bi,bj)%c(i,1), size(front%r(bi,bj)%c,1), &
          & x%front_rhs(node)%c((bj-1)*front%nb+1,1), front%n)
  end do



  ! 2nd part, multiply by the blocks in npiv+1:n
  n   = size(x%p,2)
  npj = (front%npiv)/front%nb + 1
  npi = (front%npiv-1)/front%mb + 1
  do bj=npj, front%nc
     do bi=1, npi
        k = min(bi*front%mb,front%npiv) - (bi-1)*front%mb
        j = max(1,front%npiv-(bj-1)*front%nb+1)

        m = size(front%r(bi,bj)%c,2)-j+1
        if(m.gt.0) call _xgemm(tr, 'n', &
             & m, n, k, &
             & _qrm_mone, &
             & front%r(bi,bj)%c(1,j), size(front%r(bi,bj)%c,1), &
             & x%front_rhs(node)%c((bi-1)*front%mb+1,1), front%n, &
             & _qrm_one, &
             & x%front_rhs(node)%c((bj-1)*front%nb+j,1), front%n)
     end do
  end do

  x%p(front%rows(1:front%npiv),:) = x%front_rhs(node)%c(1:front%npiv,:)
  if(front%m.gt.front%ne) x%p(front%rows(front%ne+1:front%m),:) = _qrm_zero

10 continue

9999 continue ! error
  return
end subroutine _qrm_front_rt


subroutine _qrm_assemble_rt(qrm_mat, front, b, x, info)
  use _qrm_spmat_mod
  use qrm_mem_mod
  use qrm_error_mod
  use qrm_common_mod
  use _qrm_sdata_mod
  implicit none

  type(_qrm_spmat_type), target  :: qrm_mat
  type(_qrm_front_type)          :: front
  type(_qrm_rhs_type)            :: b, x
  integer, optional              :: info

  type(qrm_adata_type), pointer  :: adata
  type(_qrm_front_type), pointer :: child
  integer                        :: n, c, node

  ! error management
  integer                        :: err
  character(len=*), parameter    :: name='qrm_assemble_rt'

  err = 0

  adata => qrm_mat%adata

  node = front%num
  n = size(b%p,2)

  if(front%ne .gt. 0) then
     call qrm_alloc(x%front_rhs(node)%c, front%n, n, err)
     __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)

     x%front_rhs(node)%c(1:front%npiv,:) = b%p(front%cols(1:front%npiv),:)
     x%front_rhs(node)%c(front%npiv+1:front%n,:) = _qrm_zero
  end if

  do c = adata%childptr(node), adata%childptr(node+1)-1
     child => qrm_mat%fdata%front_list(adata%child(c))

     ! the double condition below should be redundant but you never know...
     if((child%n.gt.child%npiv) .or. (front%ne.gt.0)) &
          & x%front_rhs(node)%c(child%colmap,:) = &
          & x%front_rhs(node)%c(child%colmap,:) + &
          & x%front_rhs(child%num)%c(child%npiv+1:child%n,:)

     call qrm_dealloc(x%front_rhs(child%num)%c)
  end do

9999 continue

  if(present(info)) info = err
  return

end subroutine _qrm_assemble_rt





subroutine _qrm_apply_subtree(transp, qrm_mat, iroot, b, work, info)
  use _qrm_spmat_mod
  use qrm_common_mod
  use qrm_string_mod
  use _qrm_solve_mod, protect => _qrm_apply_subtree
  use _qrm_sdata_mod
  implicit none
  
  type(_qrm_spmat_type), intent(in), target :: qrm_mat
  type(_qrm_rhs_type)                       :: b
  _qrm_data, intent(inout)                  :: work(:,:)
  character(len=*), intent(in)              :: transp
  integer                                   :: iroot
  integer, optional                         :: info

  integer                                   :: inode, node, f, root, ileaf
  type(qrm_adata_type), pointer             :: adata
  type(_qrm_fdata_type), pointer            :: fdata
  type(_qrm_front_type), pointer            :: front

  integer                           :: err
  character(len=*), parameter       :: name='qrm_apply_subtree'

  err = 0
  
  ! simplify
  adata => qrm_mat%adata
  fdata => qrm_mat%fdata
  root  = adata%torder(iroot) 

  if(qrm_str_tolower(transp(1:1)) .eq. _qrm_transp) then
     inode = adata%small(root)
     subtreeqt: do

        node = adata%torder(inode)
        front => fdata%front_list(node)
        call _qrm_assemble_qt(qrm_mat, front, b, err)
        __QRM_INFO_CHECK(err, name,'qrm_assemble_qt',9999)
        call _qrm_front_qt(front, b, work)
        
        if(front%num.eq.root) exit subtreeqt
        inode = inode+1
        
     end do subtreeqt
  else
     inode = iroot
     ileaf = adata%small(root)
     
     subtreeq: do

        node = adata%torder(inode)
        front => fdata%front_list(node)
        call _qrm_front_q(front, b, work)
        call _qrm_assemble_q(qrm_mat, front, b, err)
        __QRM_INFO_CHECK(err, name,'qrm_assemble_q',9999)
        
        if(inode.eq.ileaf) exit subtreeq
        inode = inode-1
        
     end do subtreeq
     
  end if

9999 continue

  if(present(info)) info = err
  return

end subroutine _qrm_apply_subtree


subroutine _qrm_solve_subtree(transp, qrm_mat, iroot, b, x, info)
  use _qrm_spmat_mod
  use qrm_common_mod
  use qrm_string_mod
  use _qrm_solve_mod, protect => _qrm_solve_subtree
  use _qrm_sdata_mod
  implicit none
  
  type(_qrm_spmat_type), intent(in), target :: qrm_mat
  type(_qrm_rhs_type)                       :: b, x
  character(len=*), intent(in)              :: transp
  integer                                   :: iroot
  integer, optional                         :: info

  integer                                   :: inode, node, f, root, ileaf
  type(qrm_adata_type), pointer             :: adata
  type(_qrm_fdata_type), pointer            :: fdata
  type(_qrm_front_type), pointer            :: front

  integer                           :: err
  character(len=*), parameter       :: name='qrm_apply_subtree'

  err = 0

  ! simplify
  adata => qrm_mat%adata
  fdata => qrm_mat%fdata
  root  = adata%torder(iroot) 

  if(qrm_str_tolower(transp(1:1)) .eq. _qrm_transp) then
     inode = adata%small(root)
     subtreert: do

        node = adata%torder(inode)
        front => fdata%front_list(node)
        call _qrm_assemble_rt(qrm_mat, front, b, x, err)
        __QRM_INFO_CHECK(err, name,'qrm_assemble_rt',9999)
        call _qrm_front_rt(front, b, x)
        
        if(front%num.eq.root) exit subtreert
        inode = inode+1
        
     end do subtreert
  else
     inode = iroot
     ileaf = adata%small(root)
     
     subtreer: do

        node = adata%torder(inode)
        front => fdata%front_list(node)
        call _qrm_front_r(front, b, x)
        call _qrm_assemble_r(qrm_mat, front, b, x, err)
        __QRM_INFO_CHECK(err, name,'qrm_assemble_r',9999)
        
        if(inode.eq.ileaf) exit subtreer
        inode = inode-1
        
     end do subtreer
     
  end if

9999 continue

  if(present(info)) info = err
  return

end subroutine _qrm_solve_subtree
