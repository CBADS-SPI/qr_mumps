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
!> @file qrm_starpu_cpu_funcs.F90
!! This file holds the CPU functions for StarPU codelets
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


! #define prnt_tsk


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  GEQRT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine _qrm_starpu_geqrt_cpu_func(buffers, cl_arg) bind(C)
  use iso_c_binding
  use starpu_f_mod
  use _qrm_factorization_mod
  use _qrm_starpu_codelets_mod
  implicit none

  type(c_ptr), value        :: cl_arg
  type(c_ptr), value        :: buffers

  integer, pointer          :: stair(:)
  type(c_ptr), target       :: cstair, a, t, w
  integer, target           :: bk, bi, mbs, nbs, ib, lda, ldt, ma, na, mt, nt, i, j, mw, nw, ldw
  _qrm_data, pointer        :: aik(:,:), tik(:,:), work(:,:)
  integer                   :: m ,n, ofs, ofsa, ne, info, id


  call starpu_f_get_buffer(buffers, 0, c_loc(a), c_loc(ma), c_loc(na), c_loc(lda))
  call c_f_pointer(a, aik,(/lda,na/))

  call starpu_f_get_buffer(buffers, 1, c_loc(t), c_loc(mt), c_loc(nt), c_loc(ldt))
  call c_f_pointer(t, tik,(/ldt,nt/))

  call starpu_f_get_buffer(buffers, 2, c_loc(w), c_loc(mw), c_loc(nw), c_loc(ldw))
  call c_f_pointer(w, work,(/ldw,nw/))

  call _qrm_starpu_unpack_args_geqrt(cl_arg, &
       & c_loc(bk), &
       & c_loc(bi), &
       & c_loc(mbs), &
       & c_loc(nbs), &
       & c_loc(ib), &
       & c_loc(cstair))

  call c_f_pointer(cstair, stair,(/na/))

  ofs = max((bk-1)*nbs +1 - (bi-1)*mbs,1)
  ofsa = (bi-1)*mbs + ofs-1
  id = starpu_f_worker_get_id()

  m = ma - ofs+1
  n = na
  ne = min(m,n)
  ib = min(ib,ne)

  call _qrm_geqrt(m, n, ib, &
       & stair, ofsa,       &
       & aik(ofs,1), lda,   &
       & tik, ldt,          &
       & work, info)

  ! For those blocks where _geqrt is done, after the
  ! _geqrt operation the v vectors are copied into the
  ! front%v block along with the T matrices. This trick
  ! allows us to remove a false dependency between _gemqrt
  ! and _tpqrt tasks. Padding is added to avoid TRSMs on
  ! the GPU and thus the T matrices are literally stored
  ! on top of the V blocks:
  !   |TTTT|
  !   |V000|
  !   |VV00|
  !   |VVV0|
  !   |VVVV|
  i = 1
  do j=1,na
     if(i.gt.ne) exit
     if(stair(j)-ofsa .le. 0) cycle
     tik(ib+i+1:ib+m,j)=aik(ofs+i:ofs+m-1,j)
     i = i+1
  end do

#if defined(prnt_tsk)
  write(*,'(i2," GEQRT  : ",i4,2x,i4)')id,bk,bi
#endif
  return

end subroutine _qrm_starpu_geqrt_cpu_func



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  GEMQRT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine _qrm_starpu_gemqrt_cpu_func(buffers, cl_arg) bind(C)
  use iso_c_binding
  use starpu_f_mod
  use _qrm_factorization_mod
  use _qrm_starpu_codelets_mod
  implicit none

  type(c_ptr), value :: cl_arg
  type(c_ptr), value :: buffers

  integer, pointer          :: stair(:)
  type(c_ptr), target       :: cstair, cflops, v, t, c, w
  integer, target           :: bk, bi, bj, mbs, nbs, ib
  integer, target           :: ldv, ldt, ldc, ldw
  integer, target           :: mv, nv, mt, nt, mc, nc, mw, nw
  _qrm_data, pointer        :: vik(:,:), tik(:,:), cij(:,:), work(:,:)
  integer                   :: m ,n, k, ofs, ofsv, ne, info, id

  call starpu_f_get_buffer(buffers, 0, c_loc(t), c_loc(mt), c_loc(nt), c_loc(ldt))
  call c_f_pointer(t, tik,(/ldt,nt/))

  call starpu_f_get_buffer(buffers, 1, c_loc(c), c_loc(mc), c_loc(nc), c_loc(ldc))
  call c_f_pointer(c, cij,(/ldc,nc/))

  call starpu_f_get_buffer(buffers, 2, c_loc(w), c_loc(mw), c_loc(nw), c_loc(ldw))
  call c_f_pointer(w, work,(/ldw,nw/))

  call _qrm_starpu_unpack_args_gemqrt(cl_arg, &
       & c_loc(bk),                              &
       & c_loc(bi),                              &
       & c_loc(bj),                              &
       & c_loc(mbs),                             &
       & c_loc(nbs),                             &
       & c_loc(ib),                              &
       & c_loc(cstair))

  call c_f_pointer(cstair, stair,(/nt/))

  ofs = max((bk-1)*nbs +1 - (bi-1)*mbs,1)
  ofsv = (bi-1)*mbs + ofs-1
  id = starpu_f_worker_get_id()

  m = mt - ib ! front%t(bk,bj)%c has been allocated to be of the good size
  n = nc
  k = nt
  ne = min(m,k)
  ib = min(ib,ne)

  call _qrm_gemqrt('l', _qrm_transp,    &
       & m, n, k, ib,           &
       & stair, ofsv,           &
       & tik(ib+1,1), ldt,      &
       & tik(1,1), ldt,         &
       & cij(ofs,1), ldc,       &
       & work, info)

#if defined(prnt_tsk)
  write(*,'(i2," GEMQRT : ",i4,2x,i4,2x,i4)')id,bk,bi,bj
#endif
  return

end subroutine _qrm_starpu_gemqrt_cpu_func



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  TPQRT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine _qrm_starpu_tpqrt_cpu_func(buffers, cl_arg) bind(C)
  use iso_c_binding
  use starpu_f_mod
  use _qrm_factorization_mod
  use _qrm_starpu_codelets_mod
  implicit none

  type(c_ptr), value             :: cl_arg
  type(c_ptr), value             :: buffers

  integer, pointer               :: stair(:)
  type(c_ptr), target            :: cstair, a, t, b, w
  integer, target                :: bk, bi, bl, mbs, nbs, ib
  integer, target                :: ma, na, lda, mt, nt, ldt
  integer, target                :: mb, nb, ldb, mw, nw, ldw
  character(kind=c_char), target :: ts
  _qrm_data, pointer             :: aik(:,:), tlk(:,:), blk(:,:), work(:,:)
  integer                        :: m ,n, l, ofs, ofsa, ofsb, ne, info, id

  call starpu_f_get_buffer(buffers, 0, c_loc(a), c_loc(ma), c_loc(na), c_loc(lda))
  call c_f_pointer(a, aik,(/lda,na/))

  call starpu_f_get_buffer(buffers, 1, c_loc(b), c_loc(mb), c_loc(nb), c_loc(ldb))
  call c_f_pointer(b, blk,(/ldb,nb/))

  call starpu_f_get_buffer(buffers, 2, c_loc(t), c_loc(mt), c_loc(nt), c_loc(ldt))
  call c_f_pointer(t, tlk,(/ldt,nt/))

  call starpu_f_get_buffer(buffers, 3, c_loc(w), c_loc(mw), c_loc(nw), c_loc(ldw))
  call c_f_pointer(w, work,(/ldw,nw/))

  call _qrm_starpu_unpack_args_tpqrt(cl_arg, &
       & c_loc(bk),     &
       & c_loc(bi),     &
       & c_loc(bl),     &
       & c_loc(mbs),    &
       & c_loc(nbs),    &
       & c_loc(ib),     &
       & c_loc(ts),     &
       & c_loc(cstair))

  call c_f_pointer(cstair, stair,(/na/))

  ofs = max((bk-1)*nbs +1 - (bi-1)*mbs,1)
  ofsa = (bi-1)*mbs + ofs-1
  ofsb = (bl-1)*mbs

  id = starpu_f_worker_get_id()

  if(ts.eq.'t') then
     m = min(mb,nb)
     l = m
  else if(ts.eq.'s') then
     m = mb
     l = 0
  end if

  n = na
  ib = min(ib,n)

  call _qrm_tpqrt( m, n, l, ib,   &
       & stair, ofsa, ofsb,       &
       & aik(ofs,1), lda,         &
       & blk, ldb,                &
       & tlk, ldt,                &
       & work, info)

#if defined(prnt_tsk)
  write(*,'(i2," TPQRT  : ",i4,2x,i4,2x,i4)')id,bk,bi,bl
#endif
  return

end subroutine _qrm_starpu_tpqrt_cpu_func


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  TPMQRT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine _qrm_starpu_tpmqrt_cpu_func(buffers, cl_arg) bind(C)
  use iso_c_binding
  use starpu_f_mod
  use _qrm_factorization_mod
  use _qrm_starpu_codelets_mod
  implicit none

  type(c_ptr), value :: cl_arg
  type(c_ptr), value :: buffers

  integer, pointer               :: stair(:)
  type(c_ptr), target            :: cstair, v, a, t, b, w
  integer, target                :: bk, bi, bl, bj, mbs, nbs, ib
  character(kind=c_char), target :: ts
  integer, target                :: ma, na, lda, mt, nt, ldt, mb, nb, ldb
  integer, target                :: mv, nv, ldv, mw, nw, ldw
  _qrm_data, pointer             :: vlk(:,:), tlk(:,:), aij(:,:), blj(:,:), work(:,:)
  integer                        :: m ,n, k, l, ofs, ofsa, ofsb, ne, info, id

  call starpu_f_get_buffer(buffers, 0, c_loc(v), c_loc(mv), c_loc(nv), c_loc(ldv))
  call c_f_pointer(v, vlk,(/ldv,nv/))

  call starpu_f_get_buffer(buffers, 1, c_loc(t), c_loc(mt), c_loc(nt), c_loc(ldt))
  call c_f_pointer(t, tlk,(/ldt,nt/))

  call starpu_f_get_buffer(buffers, 2, c_loc(a), c_loc(ma), c_loc(na), c_loc(lda))
  call c_f_pointer(a, aij,(/lda,na/))

  call starpu_f_get_buffer(buffers, 3, c_loc(b), c_loc(mb), c_loc(nb), c_loc(ldb))
  call c_f_pointer(b, blj,(/ldb,nb/))

  call starpu_f_get_buffer(buffers, 4, c_loc(w), c_loc(mw), c_loc(nw), c_loc(ldw))
  call c_f_pointer(w, work,(/ldw,nw/))

  call _qrm_starpu_unpack_args_tpmqrt(cl_arg, &
       & c_loc(bk),  &
       & c_loc(bi),  &
       & c_loc(bl),  &
       & c_loc(bj),  &
       & c_loc(mbs), &
       & c_loc(nbs), &
       & c_loc(ib),  &
       & c_loc(ts),  &
       & c_loc(cstair))

  call c_f_pointer(cstair, stair,(/nv/))

  ofs = max((bk-1)*nbs +1 - (bi-1)*mbs,1)
  ofsa = (bi-1)*mbs + ofs-1
  ofsb = (bl-1)*mbs

  id = starpu_f_worker_get_id()

  if(ts.eq.'t') then
     m = min(mv,nv)
     l = m
  else if(ts.eq.'s') then
     m = mv
     l = 0
  end if
  n = na
  k = nv
  ib = min(ib,k)

  call _qrm_tpmqrt( 'l', _qrm_transp,    &
       & m, n, k, l, ib,         &
       & stair, ofsa, ofsb,      &
       & vlk, ldv,               &
       & tlk, ldt,               &
       & aij(ofs,1), lda,        &
       & blj, ldb,               &
       & work, info)

#if defined(prnt_tsk)
  write(*,'(i2," TPMQRT : ",i4,2x,i4,2x,i4,2x,i4)')id,bk,bi,bl,bj
#endif
  return

end subroutine _qrm_starpu_tpmqrt_cpu_func


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  do_subtree
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine _qrm_starpu_do_subtree_cpu_func(buffers, cl_arg) bind(C)
  use _qrm_spmat_mod
  use _qrm_fdata_mod
  use _qrm_factorization_mod
  use _qrm_starpu_codelets_mod
  use iso_c_binding
  use qrm_dscr_mod
  implicit none

  type(c_ptr), value :: cl_arg
  type(c_ptr), value :: buffers

  type(c_ptr), target            :: qrm_mat_c, qrm_dscr_c, cflops
  integer, target                :: rootnum
  real(kind(1.d0)), target       :: flops

  type(qrm_dscr_type), pointer   :: qrm_dscr
  type(_qrm_spmat_type), pointer :: qrm_mat
  type(_qrm_front_type), pointer :: root
  integer                        :: info

  call _qrm_starpu_unpack_args_subtree(cl_arg, &
       & c_loc(qrm_dscr_c), &
       & c_loc(qrm_mat_c), &
       & c_loc(rootnum))

  call c_f_pointer(qrm_dscr_c, qrm_dscr)
  call c_f_pointer(qrm_mat_c, qrm_mat)

  ! until StarPU has a proper error handling method we have to check
  ! manually
  if(qrm_dscr%err_status.ne.0) return
  
  root => qrm_mat%fdata%front_list(rootnum)

  call _qrm_do_subtree(qrm_mat, root, flops, info)
  if(info.ne.0) qrm_dscr%err_status=info
  
  return
end subroutine _qrm_starpu_do_subtree_cpu_func



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  INIT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine _qrm_starpu_init_cpu_func(buffers, cl_arg) bind(C)
  use _qrm_spmat_mod
  use _qrm_fdata_mod
  use _qrm_factorization_mod
  use _qrm_starpu_codelets_mod
  use iso_c_binding
  use qrm_dscr_mod
  implicit none

  type(c_ptr), value             :: cl_arg
  type(c_ptr), value             :: buffers

  type(c_ptr), target            :: qrm_mat_c, qrm_dscr_c
  integer, target                :: fnum

  type(qrm_dscr_type), pointer   :: qrm_dscr
  type(_qrm_spmat_type), pointer :: qrm_mat
  type(_qrm_front_type), pointer :: front
  integer                        :: info


  call _qrm_starpu_unpack_args_init(cl_arg, &
       & c_loc(qrm_dscr_c), &
       & c_loc(qrm_mat_c), &
       & c_loc(fnum))


  call c_f_pointer(qrm_dscr_c, qrm_dscr)
  call c_f_pointer(qrm_mat_c, qrm_mat)

  ! until StarPU has a proper error handling method we have to check
  ! manually
  if(qrm_dscr%err_status.ne.0) return

  front => qrm_mat%fdata%front_list(fnum)

  call _qrm_init_front(qrm_mat, front, info)
  if(info.ne.0) qrm_dscr%err_status=info

#if defined(prnt_tsk)
  write(*,'("Init  --  inode: ",i4)')fnum
#endif
  
  return
end subroutine _qrm_starpu_init_cpu_func



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  ASSEMBLE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine _qrm_starpu_assemble_cpu_func(buffers, cl_arg) bind(C)
  use _qrm_spmat_mod
  use _qrm_fdata_mod
  use _qrm_factorization_mod
  use _qrm_starpu_codelets_mod
  use iso_c_binding
  implicit none

  type(c_ptr), value :: cl_arg
  type(c_ptr), value :: buffers

  type(c_ptr), target            :: qrm_mat_c
  integer, target                :: fnum, br, bc

  type(_qrm_spmat_type), pointer :: qrm_mat
  type(_qrm_front_type), pointer :: front
  integer                        :: err

  call _qrm_starpu_unpack_args_asm(cl_arg, &
       & c_loc(qrm_mat_c), &
       & c_loc(fnum), &
       & c_loc(br),&
       & c_loc(bc))

  call c_f_pointer(qrm_mat_c, qrm_mat)

  call _qrm_assemble_block(qrm_mat, fnum, br, bc)
#if defined(prnt_tsk)
  write(*,'("Assmb --  inode: ",i4,2x,i4,2x,i4)')fnum,br,bc
#endif
  return
end subroutine _qrm_starpu_assemble_cpu_func






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  CLEAN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine _qrm_starpu_clean_cpu_func(buffers, cl_arg) bind(C)
  use _qrm_spmat_mod
  use _qrm_fdata_mod
  use _qrm_factorization_mod
  use _qrm_starpu_codelets_mod
  use iso_c_binding
  use qrm_dscr_mod
  implicit none

  type(c_ptr), value             :: cl_arg
  type(c_ptr), value             :: buffers

  type(c_ptr), target            :: qrm_mat_c, qrm_dscr_c
  integer, target                :: fnum

  type(qrm_dscr_type), pointer   :: qrm_dscr
  type(_qrm_spmat_type), pointer :: qrm_mat
  type(_qrm_front_type), pointer :: front
  integer                        :: info

  call _qrm_starpu_unpack_args_clean(cl_arg, &
       & c_loc(qrm_dscr_c),                  &
       & c_loc(qrm_mat_c),                   &
       & c_loc(fnum))

  call c_f_pointer(qrm_dscr_c, qrm_dscr)
  call c_f_pointer(qrm_mat_c, qrm_mat)

  ! until StarPU has a proper error handling method we have to check
  ! manually
  if(qrm_dscr%err_status.ne.0) return

  front => qrm_mat%fdata%front_list(fnum)

  call _qrm_clean_front(qrm_mat, front, info)
  if(info.ne.0) qrm_dscr%err_status=info

#if defined(prnt_tsk)
  write(*,'("Clean --  inode: ",i4)')fnum
#endif
  return
end subroutine _qrm_starpu_clean_cpu_func




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  APPLY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine _qrm_starpu_apply_cpu_func(buffers, cl_arg) bind(C)
  use _qrm_spmat_mod
  use _qrm_fdata_mod
  use _qrm_sdata_mod
  use _qrm_factorization_mod
  use _qrm_starpu_codelets_mod
  use _qrm_solve_mod  ! TODO: add dep in makefile
  use iso_c_binding
  use qrm_dscr_mod
  implicit none

  type(c_ptr), value             :: cl_arg
  type(c_ptr), value             :: buffers

  type(c_ptr), target            :: qrm_mat_c, b_c, qrm_dscr_c
  integer, target                :: inode, node
  type(_qrm_front_type), pointer :: front
  type(_qrm_rhs_type), pointer   :: b

  type(qrm_dscr_type), pointer   :: qrm_dscr
  type(_qrm_spmat_type), pointer :: qrm_mat
  character(kind=c_char), target :: transp
  integer                        :: info
  type(c_ptr), target            :: w
  integer, target                :: mw, nw, ldw
  _qrm_data, pointer             :: work(:,:)
  logical                        :: issmall

  call starpu_f_get_buffer(buffers, 0, c_loc(w), c_loc(mw), c_loc(nw), c_loc(ldw))
  call c_f_pointer(w, work,(/ldw,nw/))

  call _qrm_starpu_unpack_args_apply(cl_arg, &
       & c_loc(qrm_dscr_c), &
       & c_loc(qrm_mat_c), &
       & c_loc(b_c), &
       & c_loc(inode),&
       & c_loc(transp))

  call c_f_pointer(qrm_dscr_c, qrm_dscr)
  call c_f_pointer(qrm_mat_c, qrm_mat)
  call c_f_pointer(b_c, b)

  ! until StarPU has a proper error handling method we have to check
  ! manually
  if(qrm_dscr%err_status.ne.0) return

  node = qrm_mat%adata%torder(inode)
  issmall = qrm_mat%adata%small(node) .gt. 0

  if(issmall) then
     call _qrm_apply_subtree(transp, qrm_mat, inode, b, work, info)
  else
     front => qrm_mat%fdata%front_list(node)
     if(transp .eq. _qrm_transp) then
        call _qrm_assemble_qt(qrm_mat, front, b, info)
        if(info.eq.0) call _qrm_front_qt(front, b, work)
     else
        call _qrm_front_q(front, b, work)
        call _qrm_assemble_q(qrm_mat, front, b, info)
     end if
  end if

  if(info.ne.0) qrm_dscr%err_status=info

#if defined(prnt_tsk)
  write(*,'("Apply --  inode: ",i4,"   size:",i4,2x,i4)')inode,size(b%p,1),size(b%p,2)
#endif
  return

end subroutine _qrm_starpu_apply_cpu_func


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  SOLVE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine _qrm_starpu_solve_cpu_func(buffers, cl_arg) bind(C)
  use _qrm_spmat_mod
  use _qrm_fdata_mod
  use _qrm_sdata_mod
  use _qrm_factorization_mod
  use _qrm_starpu_codelets_mod
  use _qrm_solve_mod  ! TODO: add dep in makefile
  use iso_c_binding
  use qrm_dscr_mod
  implicit none

  type(c_ptr), value             :: cl_arg
  type(c_ptr), value             :: buffers

  type(c_ptr), target            :: qrm_mat_c, qrm_dscr_c, b_c, x_c
  integer, target                :: inode, node
  type(_qrm_front_type), pointer :: front
  type(_qrm_rhs_type), pointer   :: b, x

  type(qrm_dscr_type), pointer   :: qrm_dscr
  type(_qrm_spmat_type), pointer :: qrm_mat
  character(kind=c_char), target :: transp
  integer                        :: info
  logical                        :: issmall

  call _qrm_starpu_unpack_args_solve(cl_arg, &
       & c_loc(qrm_dscr_c), &
       & c_loc(qrm_mat_c), &
       & c_loc(b_c), &
       & c_loc(x_c), &
       & c_loc(inode),&
       & c_loc(transp))

  call c_f_pointer(qrm_dscr_c, qrm_dscr)
  call c_f_pointer(qrm_mat_c, qrm_mat)
  call c_f_pointer(b_c, b)
  call c_f_pointer(x_c, x)

  ! until StarPU has a proper error handling method we have to check
  ! manually
  if(qrm_dscr%err_status.ne.0) return

  node = qrm_mat%adata%torder(inode)
  issmall = qrm_mat%adata%small(node) .gt. 0

  if(issmall) then
     call _qrm_solve_subtree(transp, qrm_mat, inode, b, x, info)
  else
     front => qrm_mat%fdata%front_list(node)
     if(transp .eq. _qrm_transp) then
        call _qrm_assemble_rt(qrm_mat, front, b, x, info)
        if(info.eq.0) call _qrm_front_rt(front, b, x)
     else
        call _qrm_front_r(front, b, x)
        call _qrm_assemble_r(qrm_mat, front, b, x, info)
     end if
  end if

  if(info.ne.0) qrm_dscr%err_status=info

#if defined(prnt_tsk)
  write(*,'("Solve --  inode: ",i4,"   size:",i4,2x,i4)')inode,size(b%p,1),size(b%p,2)
#endif
  return

end subroutine _qrm_starpu_solve_cpu_func


subroutine _qrm_starpu_analyse_cpu_func(buffers, cl_arg) bind(C)
  use _qrm_spmat_mod
  use _qrm_analysis_mod
  use iso_c_binding
  use _qrm_starpu_codelets_mod
  use qrm_dscr_mod
  implicit none

  type(c_ptr), value             :: cl_arg
  type(c_ptr), value             :: buffers

  type(c_ptr), target            :: qrm_mat_c, qrm_dscr_c
  type(_qrm_spmat_type), pointer :: qrm_mat
  type(qrm_dscr_type), pointer   :: qrm_dscr
  character(kind=c_char), target :: transp
  integer                        :: info

  call _qrm_starpu_unpack_args_analyse(cl_arg, &
       & c_loc(qrm_dscr_c), &
       & c_loc(qrm_mat_c), &
       & c_loc(transp))

  call c_f_pointer(qrm_dscr_c, qrm_dscr)
  call c_f_pointer(qrm_mat_c, qrm_mat)

  call _qrm_analysis_core(qrm_mat, qrm_dscr, transp, info)
  if(info.ne.0) qrm_dscr%err_status=info
  return

end subroutine _qrm_starpu_analyse_cpu_func
