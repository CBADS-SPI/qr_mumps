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
!> @file qrm_factorization_tasks.F90
!! This file holds the implementation of various factorization tasks
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  GEQRT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! submit geqrt task to StarPU
subroutine _qrm_geqrt_task(qrm_dscr, front, bk, bi, work, prio)
  use _qrm_spmat_mod
  use _qrm_fdata_mod
  use iso_c_binding
  use qrm_common_mod
  use _qrm_utils_mod
  use qrm_dscr_mod
  use _qrm_factorization_mod, protect => _qrm_geqrt_task
#if defined(have_starpu)
  use _qrm_starpu_codelets_mod
#endif
  implicit none

  type(qrm_dscr_type)                 :: qrm_dscr
  type(_qrm_front_type), target       :: front
  integer                             :: fnum, bk, bi
  type(_qrm_bc_type)                  :: work
  integer, optional                   :: prio

  type(c_ptr)                         :: stair_c, sym_handle
  real(kind=c_double)                 :: flops_c
  logical                             :: sym
  
  integer                             :: iprio
  integer                             :: m, n, ofs, ofsa, ne, ib, info
  integer                             :: ma, na, lda, ldt, i, j

  if( (front%n .le. 0) .or. (front%m .le. 0)) return

#if defined (have_starpu)
  stair_c   = c_loc(front%stair((bk-1)*front%nb+1))

  if(present(prio)) then
     iprio = prio
  else
     iprio = 1
  end if

  ! all the geqrt tasks on blocks that have nothing on their left
  ! should take the front symbolic handle to respect dependencies wrt
  ! the init routine.
  sym_handle = c_null_ptr
  if(bk.eq.1) then
     sym_handle = front%sym_handle
  ! 14-06-2016: not sur sure why this else was initially coded this
  ! way. maybe it's just a bug but I'm leaving it for future reference
  ! else if((bk.gt.1) .and. (bk.lt.front%nc)) then
  else
     if(((front%stair((bk-1)*front%nb)-1)/front%mb+1).lt.bi) then
        sym_handle = front%sym_handle
     end if
  end if

  if(    (front%ne.eq.front%npiv) .and. &
       & (bk.eq.front%nc)         .and. &
       & (bi.eq.front%nr)         ) then
     sym_handle=front%sym_handle
  end if
  
  call _qrm_insert_geqrt_task_c(                      &
       & qrm_dscr,                                    &
       & front%bc(bi,bk)%hdl, front%t(bi,bk)%hdl,     &
       & work%hdl,                                    &
       & sym_handle,                                  &
       & bk, bi,                                      &
       & front%mb, front%nb, front%ib,                &
       & stair_c, flops_c, iprio)
  
#else

  ma  = size(front%bc(bi,bk)%c,1)
  na  = size(front%bc(bi,bk)%c,2)
  lda = size(front%bc(bi,bk)%c,1)

  ldt = size(front%t(bi,bk)%c,1)
  
  ofs = max((bk-1)*front%nb +1 - (bi-1)*front%mb,1)
  ofsa = (bi-1)*front%mb + ofs-1
  m = ma - ofs+1
  n = na
  ne = min(m,n)
  ib = min(front%ib,ne)
  
  call _qrm_geqrt(m, n, ib,                    &
       & front%stair((bk-1)*front%nb+1), ofsa, &
       & front%bc(bi,bk)%c(ofs,1), lda,        & 
       & front%t(bi,bk)%c(1,1), ldt,           &
       & work%c(1,1), info)

  i = 1
  do j=1,na
     if(i.gt.ne) exit
     if(front%stair((bk-1)*front%nb+j)-ofsa .le. 0) cycle
     front%t(bi,bk)%c(ib+i+1:ib+m,j)=front%bc(bi,bk)%c(ofs+i:ofs+m-1,j)
     i = i+1
  end do
  
#endif

  
  return

end subroutine _qrm_geqrt_task


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  GEMQRT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine _qrm_gemqrt_task(qrm_dscr, front, bk, bi, bj, work, prio)
  use _qrm_spmat_mod
  use _qrm_fdata_mod
  use iso_c_binding
  use qrm_common_mod
  use _qrm_utils_mod
  use qrm_dscr_mod
  use _qrm_factorization_mod, protect => _qrm_gemqrt_task
#if defined(have_starpu)
  use _qrm_starpu_codelets_mod
  ! use starpu_f_mod
#endif
  implicit none

  type(qrm_dscr_type)             :: qrm_dscr
  type(_qrm_front_type), target   :: front
  integer                         :: fnum, bk, bi, bj
  type(_qrm_bc_type)              :: work
  integer, optional               :: prio

  type(c_ptr)                     :: stair_c, sym_handle
  real(kind=c_double)             :: flops_c
  real(c_double)                  :: flops_task
  integer                         :: iprio
  integer                         :: m, n, k, ofs, ofsv, ib, ne, info
  integer                         :: mt, nt, ldt, mc, nc, ldc
  
  if( (front%n .le. 0) .or. (front%m .le. 0)) return
  

#if defined (have_starpu)
  stair_c   = c_loc(front%stair((bk-1)*front%nb+1))
  
  if(present(prio)) then
     iprio = prio
  else
     iprio = 1
  end if

  sym_handle = c_null_ptr
  if(    (front%ne.eq.front%npiv) .and. & ! only if no asm
       & (bk.eq.front%np)         .and. &
       & (bj.gt.front%np)         .and. &
       & (bi.eq.front%nr)         ) then
     sym_handle=front%sym_handle
  end if
  
  call _qrm_insert_gemqrt_task_c(        &
       & qrm_dscr,                       &
       & front%t(bi,bk)%hdl,             &
       & front%bc(bi,bj)%hdl,            &
       & work%hdl,                       &
       & sym_handle,                     &
       & bk, bi, bj,                     &
       & front%mb, front%nb, front%ib,   &
       & stair_c, flops_c, iprio)
#else
  
  mc  = size(front%bc(bi,bj)%c,1)
  nc  = size(front%bc(bi,bj)%c,2)
  ldc = size(front%bc(bi,bj)%c,1)
  
  mt  = size(front%t(bi,bk)%c,1)
  nt  = size(front%t(bi,bk)%c,2)
  ldt = size(front%t(bi,bk)%c,1)
  
  ofs = max((bk-1)*front%nb +1 - (bi-1)*front%mb,1)
  ofsv = (bi-1)*front%mb + ofs-1

  m = mt - front%ib ! front%t(bk,bj)%c has been allocated to be of the good size
  n = nc
  k = nt
  ne = min(m,k)
  ib = min(front%ib,ne)
  ! starpu not available, execute task

  call _qrm_gemqrt('l', _qrm_transp,                                  &
       & m, n, k, ib,                                         & 
       & front%stair((bk-1)*front%nb+1), ofsv,                &
       & front%t(bi,bk)%c(ib+1,1), ldt , &
       & front%t(bi,bk)%c(1,1),    ldt , &
       & front%bc(bi,bj)%c(ofs,1), ldc , &
       & work%c(1,1), info)

#endif

  return

end subroutine _qrm_gemqrt_task

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  TPQRT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine _qrm_tpqrt_task(qrm_dscr, front, bk, bi, bl, ts, work, prio)
  use _qrm_spmat_mod
  use _qrm_fdata_mod
  use iso_c_binding
  use qrm_common_mod
  use _qrm_utils_mod
  use qrm_dscr_mod
  use _qrm_factorization_mod, protect => _qrm_tpqrt_task
#if defined(have_starpu)
  use _qrm_starpu_codelets_mod
  ! use starpu_f_mod
#endif
  implicit none

  type(qrm_dscr_type)             :: qrm_dscr
  type(_qrm_front_type), target   :: front
  integer                         :: fnum, bk, bi, bl
  character(kind=c_char), value   :: ts
  type(_qrm_bc_type)              :: work
  integer, optional               :: prio
  
  type(c_ptr)                     :: stair_c, sym_handle
  real(kind=c_double)             :: flops_c
  integer                         :: iprio
  integer                         :: m, n, k, l, ofs, ofsa, ofsb, ib, info, ft, lt, nsd
  integer                         :: mt, nt, ldt, ma, na, lda, mb, nb, ldb

  if( (front%n .le. 0) .or. (front%m .le. 0)) return


#if defined (have_starpu)
  stair_c   = c_loc(front%stair((bk-1)*front%nb+1))

  if(present(prio)) then
     iprio = prio
  else
     iprio = 1
  end if

  sym_handle = c_null_ptr
  ft = ((bk-1)*front%nb)/front%mb+1
  lt = (front%stair(min(bk*front%nb, front%n))-1)/front%mb+1
  
  if(    (front%ne.eq.front%npiv) .and. & ! only if no asm
       & (bk.eq.front%nc)         .and. &
       & (bi.eq.ft) ) then
     
     nsd = (lt-ft)/front%bh+1
     if((nsd.eq.1) .and. (bl.eq.lt)) then
        sym_handle=front%sym_handle
     else if((nsd.gt.1) .and. &
          &  (ft+(bl-bi)*2.gt.lt )) then
        sym_handle=front%sym_handle
     end if
  end if
  
  call _qrm_insert_tpqrt_task_c(                       &
       & qrm_dscr,                                     &
       & front%bc(bi,bk)%hdl, front%bc(bl,bk)%hdl,     &
       & front%t2(bl,bk)%hdl,                          &
       & work%hdl,                                     &
       & sym_handle,                                   &
       & bk, bi, bl,                                   &
       & front%mb, front%nb, front%ib,                 &
       & ts, stair_c, flops_c, iprio)
#else

  ma  = size(front%bc(bi,bk)%c,1)
  na  = size(front%bc(bi,bk)%c,2)
  lda = size(front%bc(bi,bk)%c,1)
  
  mb  = size(front%bc(bl,bk)%c,1)
  nb  = size(front%bc(bl,bk)%c,2)
  ldb = size(front%bc(bl,bk)%c,1)
  
  mt  = size(front%t2(bl,bk)%c,1)
  nt  = size(front%t2(bl,bk)%c,2)
  ldt = size(front%t2(bl,bk)%c,1)

  ofs = max((bk-1)*front%nb +1 - (bi-1)*front%mb,1)
  ofsa = (bi-1)*front%mb + ofs-1
  ofsb = (bl-1)*front%mb

  if(ts.eq.'t') then
     m = min(mb,nb)
     l = m
  else if(ts.eq.'s') then
     m = mb
     l = 0
  end if

  n = na
  ib = min(front%ib,n)
  
  call _qrm_tpqrt( m, n, l, ib,                      &
       & front%stair((bk-1)*front%nb+1), ofsa, ofsb, &
       & front%bc(bi,bk)%c(ofs,1), lda,              &
       & front%bc(bl,bk)%c(1,1),   ldb,              &
       & front%t2(bl,bk)%c(1,1),   ldt,              &
       & work%c(1,1), info)

#endif

  return

end subroutine _qrm_tpqrt_task

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  TPMQRT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine _qrm_tpmqrt_task(qrm_dscr, front, bk, bi, bl, bj, ts, work, prio)
  use _qrm_spmat_mod
  use _qrm_fdata_mod
  use iso_c_binding
  use qrm_common_mod
  use _qrm_utils_mod
  use qrm_dscr_mod
  use _qrm_factorization_mod, protect => _qrm_tpmqrt_task
#if defined(have_starpu)
  use _qrm_starpu_codelets_mod
#endif
  implicit none

  type(qrm_dscr_type)             :: qrm_dscr
  type(_qrm_front_type), target   :: front
  integer                         :: fnum, bk, bi, bl, bj
  character(kind=c_char), value   :: ts
  type(_qrm_bc_type)              :: work
  integer, optional               :: prio

  type(c_ptr)                     :: stair_c, sym_handle
  real(kind=c_double)             :: flops_c
  integer                         :: iprio
  integer                         :: m, n, k, l, ofs, ofsa, ofsb, ib, info, ft, lt, nsd
  integer                         :: mt, nt, ldt, ma, na, lda, mb, nb, ldb, mv, nv, ldv
  
  if( (front%n .le. 0) .or. (front%m .le. 0)) return

#if defined (have_starpu)
  stair_c = c_loc(front%stair((bk-1)*front%nb+1))

  if(present(prio)) then
     iprio = prio
  else
     iprio = 1
  end if

  sym_handle = c_null_ptr
  ft = ((bk-1)*front%nb)/front%mb+1
  lt = (front%stair(min(bk*front%nb, front%n))-1)/front%mb+1
  
  if(    (front%ne.eq.front%npiv) .and. & ! only if no asm
       & (bk.eq.front%np)         .and. &
       & (bj.gt.front%np)         .and. &
       & (bi.eq.ft) ) then
     
     nsd = (lt-ft)/front%bh+1
     if((nsd.eq.1) .and. (bl.eq.lt)) then
        sym_handle=front%sym_handle
     else if((nsd.gt.1) .and. &
          &  (ft+(bl-bi)*2.gt.lt )) then
        sym_handle=front%sym_handle
     end if
  end if

  call _qrm_insert_tpmqrt_task_c(                      &
       & qrm_dscr,                                     &
       & front%bc(bl,bk)%hdl, front%t2(bl,bk)%hdl,     &
       & front%bc(bi,bj)%hdl, front%bc(bl,bj)%hdl,     &
       & work%hdl,                                     &
       & sym_handle,                                   &
       & bk, bi, bl, bj,                               &
       & front%mb, front%nb, front%ib,                 &
       & ts, stair_c, flops_c, iprio)
#else

  ma  = size(front%bc(bi,bj)%c,1)
  na  = size(front%bc(bi,bj)%c,2)
  lda = size(front%bc(bi,bj)%c,1)
  
  mb  = size(front%bc(bl,bj)%c,1)
  nb  = size(front%bc(bl,bj)%c,2)
  ldb = size(front%bc(bl,bj)%c,1)
  
  mv  = size(front%bc(bl,bk)%c,1)
  nv  = size(front%bc(bl,bk)%c,2)
  ldv = size(front%bc(bl,bk)%c,1)
  
  mt  = size(front%t2(bl,bk)%c,1)
  nt  = size(front%t2(bl,bk)%c,2)
  ldt = size(front%t2(bl,bk)%c,1)

  ofs = max((bk-1)*front%nb +1 - (bi-1)*front%mb,1)
  ofsa = (bi-1)*front%mb + ofs-1
  ofsb = (bl-1)*front%mb

  if(ts.eq.'t') then
     m = min(mv,nv)
     l = m
  else if(ts.eq.'s') then
     m = mv
     l = 0
  end if
  n = na
  k = nv
  ib = min(front%ib,k)
  
  call _qrm_tpmqrt( 'l', _qrm_transp,                        &
       & m, n, k, l, ib,                             &
       & front%stair((bk-1)*front%nb+1), ofsa, ofsb, &
       & front%bc(bl,bk)%c(1,1),   ldv,              &
       & front%t2(bl,bk)%c(1,1),   ldt,              &
       & front%bc(bi,bj)%c(ofs,1), lda,              &
       & front%bc(bl,bj)%c(1,1),   ldb,              &
       & work%c(1,1), info)

#endif
  
  return

end subroutine _qrm_tpmqrt_task



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  do_subtree
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine _qrm_do_subtree_task(qrm_dscr, qrm_mat, rootnum, prio, info)
  use _qrm_spmat_mod
  use _qrm_fdata_mod
  use iso_c_binding
  use qrm_dscr_mod
  use _qrm_factorization_mod, protect=>_qrm_do_subtree_task
#if defined(have_starpu)
  use _qrm_starpu_codelets_mod
#endif
  implicit none
  
  type(qrm_dscr_type)                 :: qrm_dscr
  type(_qrm_spmat_type), target       :: qrm_mat
  integer                             :: rootnum
  integer, optional                   :: prio, info
  
  type(_qrm_front_type), pointer      :: father, front
  integer                             :: f
  type(c_ptr)                         :: qrm_mat_c, father_sym_handle_c
  real(kind=c_double)                 :: flops_c
  real(kind(1.d0))                    :: rflops
  integer                             :: iprio

  ! error management
  integer                             :: err
  character(len=*), parameter         :: name='qrm_do_subtree_task'

  err = 0
  
#if defined(have_starpu)
  front => qrm_mat%fdata%front_list(rootnum)
     
  qrm_mat_c = c_loc(qrm_mat)
  if(present(prio)) then
     iprio = prio
  else
     iprio = 1
  end if

  call _qrm_insert_do_subtree_task_c(qrm_dscr, &
       & qrm_mat_c,                            &
       & int(rootnum, kind=c_int),             &
       & front%sym_handle,                     &
       & flops_c, iprio)
#else
  call _qrm_do_subtree(qrm_mat, qrm_mat%fdata%front_list(rootnum), rflops, err)
#endif

  if(present(info)) info = err
  return
  
end subroutine _qrm_do_subtree_task

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  INIT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine _qrm_init_task(qrm_dscr, qrm_mat, fnum, prio, info)
  ! the Intel compler bugs here when used with vectorization options
  !DIR$ NOOPTIMIZE  
  
  use iso_c_binding
  use _qrm_fdata_mod
  use _qrm_spmat_mod
  use qrm_dscr_mod
  use _qrm_factorization_mod, protect=>_qrm_init_task
#if defined(have_starpu)
  use _qrm_starpu_codelets_mod
#endif
  implicit none

  type(qrm_dscr_type)            :: qrm_dscr
  type(_qrm_spmat_type), target  :: qrm_mat
  integer                        :: fnum
  integer, optional              :: prio, info

  type(_qrm_front_type), pointer :: front, father
  type(c_ptr)                    :: qrm_mat_c, father_sym_handle_c
  integer                        :: f
  integer                        :: iprio
#if defined (have_starpu)
  type(c_ptr), allocatable       :: chandles(:)
  integer                        :: nc, c, i
#endif  

  ! error management
  integer                        :: err
  character(len=*), parameter    :: name='qrm_init_task'

  err = 0
  
  
#if defined (have_starpu)
  qrm_mat_c = c_loc(qrm_mat)  

  if(present(prio)) then
     iprio = prio
  else
     iprio = 1
  end if

  front => qrm_mat%fdata%front_list(fnum)

  nc = qrm_mat%adata%childptr(fnum+1)-qrm_mat%adata%childptr(fnum)
  allocate(chandles(nc))
  do i=1, nc 
     c = qrm_mat%adata%child(i-1 + qrm_mat%adata%childptr(fnum))
     chandles(i) = qrm_mat%fdata%front_list(c)%sym_handle
  end do

  call _qrm_insert_init_task_c(qrm_dscr,  &
       & qrm_mat_c,                       &
       & int(fnum, kind=c_int),           &
       & front%sym_handle,                &
       & chandles,                        &
       & nc, iprio)

  deallocate(chandles)

#else
  call _qrm_init_front(qrm_mat, qrm_mat%fdata%front_list(fnum), err)
#endif

  if(present(info)) info = err
  return

end subroutine _qrm_init_task

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  ASSEMBLE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! submit assembly block task in StarPU
subroutine _qrm_assemble_block_task(qrm_dscr, qrm_mat, fnum, br, bc, prio)
  use iso_c_binding
  use _qrm_spmat_mod
  use _qrm_fdata_mod
  use qrm_common_mod
  use qrm_dscr_mod
  use _qrm_factorization_mod, protect=>_qrm_assemble_block_task
#if defined(have_starpu)
  use _qrm_starpu_codelets_mod
#endif
  implicit none

  type(qrm_dscr_type)            :: qrm_dscr
  type(_qrm_spmat_type), target  :: qrm_mat
  integer                        :: br, bc, fnum
  integer, optional              :: prio
  
  integer                        :: f, cc, cr, fcr, fcc, fbr, fbc, lcc, i, j, nh
  type(_qrm_front_type), pointer :: front, father
  type(c_ptr)                    :: qrm_mat_c
  type(c_ptr), allocatable       :: bhandles(:)
  real(kind(1.d0))               :: t_ins_start, t_ins ! trace
  integer                        :: iprio
  
#if defined (have_starpu)
  front => qrm_mat%fdata%front_list(fnum)
  f = qrm_mat%adata%parent(front%num)
  father => qrm_mat%fdata%front_list(f)

  qrm_mat_c = c_loc(qrm_mat)
  nh = 0
  allocate(bhandles(father%nr*father%nc))

  fcr = (br-1)*front%mb
  fcc = (bc-1)*front%nb
  lcc = min(bc*front%nb,front%ne)-front%npiv

  fbc = 0
  ! The loop on j is in reverse order for handling blocks that are traversed by the diagonal
  jloop: do j=size(front%bc(br,bc)%c,2), 1, -1
  ! jloop: do j=1, size(front%bc(br,bc)%c,2)
     cc = fcc+j-front%npiv
     if(cc.le.0) cycle jloop
     if(fbc.eq.(front%colmap(cc)-1)/father%nb+1) cycle jloop
     fbc = (front%colmap(cc)-1)/father%nb+1
     fbr = 0
     
     iloop: do i=1, size(front%bc(br,bc)%c,1)
        cr = fcr+i-front%npiv
        if(cr.le.0) cycle iloop                                 ! loop if row < npiv
        if((cr.gt.lcc) .or. (fcr+i .gt. fcc+j)) cycle jloop     ! loop if below the CB
        if(fbr.eq.(front%rowmap(cr)-1)/father%mb+1) cycle iloop ! loop if within the same parent block
        fbr = (front%rowmap(cr)-1)/father%mb+1
        nh = nh+1
        bhandles(nh) = father%bc(fbr,fbc)%hdl
        ! write(*,'(i2," -- c:(",i2,",",i2,") -> f:(",i2,",",i2,")")')front%num,br,bc,fbr,fbc
     end do iloop
  end do jloop

  if(present(prio)) then
     iprio = prio
  else
     iprio = 1
  end if

  call _qrm_insert_asm_blk_task_c(qrm_dscr, &
       & qrm_mat_c, fnum,                   &
       & front%bc(br,bc)%hdl,               &  
       & front%sym_handle,                  &
       & bhandles,                          &  
       & father%sym_handle,                 &
       & nh, br, bc, iprio)

  deallocate(bhandles)
#else
  call _qrm_assemble_block(qrm_mat, fnum, br, bc)
#endif
  
  return 

end subroutine _qrm_assemble_block_task

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  CLEAN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine _qrm_clean_task(qrm_dscr, qrm_mat, fnum, prio, info)
  
  use iso_c_binding
  use _qrm_fdata_mod
  use _qrm_spmat_mod
  use qrm_dscr_mod
  use _qrm_factorization_mod, protect=>_qrm_clean_task
#if defined(have_starpu)
  use _qrm_starpu_codelets_mod
#endif
  implicit none

  type(qrm_dscr_type)            :: qrm_dscr
  type(_qrm_spmat_type), target  :: qrm_mat
  integer                        :: fnum
  integer, optional              :: prio, info

  type(_qrm_front_type), pointer :: front
  type(c_ptr)                    :: qrm_mat_c
  integer                        :: iprio

  ! error management
  integer                        :: err
  character(len=*), parameter    :: name='qrm_clean_task'

  err = 0
  
  front => qrm_mat%fdata%front_list(fnum)

#if defined (have_starpu)
  qrm_mat_c = c_loc(qrm_mat)  
  if(present(prio)) then
     iprio = prio
  else
     iprio = 1
  end if

  call _qrm_insert_clean_task_c(qrm_dscr, &
       & qrm_mat_c,                       &
       & int(fnum, kind=c_int),           &
       & front%sym_handle,                &
       & iprio)

#else
  call _qrm_clean_front(qrm_mat, front)
#endif

  if(present(info)) info = err
  return

end subroutine _qrm_clean_task

