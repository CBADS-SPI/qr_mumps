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

module _qrm_starpu_codelets_mod
  use iso_c_binding

  ! interoperable global variables
  type(c_ptr), bind(c)               :: task_duration_worker_c
  type(c_ptr), bind(c)               :: flops_worker_c
  type(c_ptr), bind(c)               :: work_c

  integer, parameter                 :: qrm_assembled_      = 0
  integer, parameter                 :: qrm_partitioned_    = 1



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! C routines interfaces
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  interface
     subroutine _qrm_insert_geqrt_task_c(     &
          & qrm_dscr,                         &
          & a_handle, t_handle,               &
          & work_handle, sym_handle,          &
          & bk, bi,                           &
          & mb, nb, ib,                       &
          & stair, flops, prio) bind(c)
       use qrm_dscr_mod
       use iso_c_binding
       type(qrm_dscr_type)           :: qrm_dscr
       type(c_ptr), value            :: a_handle, t_handle, sym_handle, work_handle
       integer(c_int), value         :: bk, bi, mb, nb, ib, prio
       type(c_ptr), value            :: stair
       real(kind=c_double), value    :: flops
     end subroutine _qrm_insert_geqrt_task_c
  end interface

  interface
     subroutine _qrm_insert_gemqrt_task_c(   &
          & qrm_dscr,                        &
          & t_handle, c_handle,              &
          & work_handle, sym_handle,         &
          & bk, bi, bj,                      &
          & mb, nb, ib,                      &
          & stair, flops, prio) bind(c)
       use iso_c_binding
       use qrm_dscr_mod
       type(qrm_dscr_type)           :: qrm_dscr
       type(c_ptr), value            :: t_handle, c_handle, work_handle, sym_handle
       integer(c_int), value         :: bk, bi, bj, mb, nb, ib, prio
       type(c_ptr), value            :: stair
       real(kind=c_double), value    :: flops
     end subroutine _qrm_insert_gemqrt_task_c
  end interface
  
  interface
     subroutine _qrm_insert_tpqrt_task_c( &
          & qrm_dscr,                     &
          & a_handle, b_handle,           &
          & t_handle, work_handle,        &
          & sym_handle,                   &
          & bk, bi, bl,                   &
          & mb, nb, ib,                   &
          & ts, stair,                    &
          & flops, prio) bind(c)
       use iso_c_binding
       use qrm_dscr_mod
       type(qrm_dscr_type)           :: qrm_dscr
       type(c_ptr), value            :: a_handle, b_handle, t_handle, work_handle, sym_handle
       integer(c_int), value         :: bk, bi, bl, mb, nb, ib, prio
       character(kind=c_char), value :: ts
       type(c_ptr), value            :: stair
       real(kind=c_double), value    :: flops
     end subroutine _qrm_insert_tpqrt_task_c
  end interface
  
  interface
     subroutine _qrm_insert_tpmqrt_task_c( &
          & qrm_dscr,                      &
          & v_handle, t_handle,            &
          & a_handle, b_handle,            &
          & work_handle, sym_handle,       &
          & bk, bi, bl, bj,                &
          & mb, nb, ib,                    &
          & ts, stair,                     &
          & flops, prio) bind(c)
       use iso_c_binding
       use qrm_dscr_mod
       type(qrm_dscr_type)           :: qrm_dscr
       type(c_ptr), value            :: a_handle, b_handle, v_handle, t_handle, work_handle, sym_handle
       integer(c_int), value         :: bk, bi, bl, bj, mb, nb, ib, prio
       character(kind=c_char), value :: ts
       type(c_ptr), value            :: stair
       real(kind=c_double), value    :: flops
     end subroutine _qrm_insert_tpmqrt_task_c
  end interface
  
  interface
     subroutine _qrm_insert_do_subtree_task_c(qrm_dscr, qrm_mat, rootnum, &
          & sym_handle, flops, prio) bind(c)
       use iso_c_binding
       use qrm_dscr_mod
       type(qrm_dscr_type)           :: qrm_dscr
       type(c_ptr), value            :: qrm_mat, sym_handle
       integer(c_int), value         :: rootnum, prio
       real(kind=c_double), value    :: flops
     end subroutine _qrm_insert_do_subtree_task_c
  end interface
  
  interface
     subroutine _qrm_insert_init_task_c(qrm_dscr, qrm_mat, fnum, sym_handle, &
          & kids_handles, nc, prio) bind(c)
       use iso_c_binding
       use qrm_dscr_mod
       type(qrm_dscr_type)           :: qrm_dscr
       type(c_ptr), value            :: qrm_mat, sym_handle
       integer(c_int), value         :: fnum, nc, prio
       type(c_ptr)                   :: kids_handles(*)
     end subroutine _qrm_insert_init_task_c
  end interface
    
  interface
     subroutine _qrm_insert_asm_blk_task_c(qrm_dscr, qrm_mat, fnum, &
          & block_handle, sym_handle,                     &
          & father_handles, father_sym_handle,            &
          & nh, br, bc, prio) bind(c)
       use iso_c_binding
       use qrm_dscr_mod
       type(qrm_dscr_type)           :: qrm_dscr
       type(c_ptr), value            :: qrm_mat, block_handle
       type(c_ptr), value            :: sym_handle, father_sym_handle
       type(c_ptr)                   :: father_handles(*)
       integer(c_int), value         :: br, bc, nh, fnum, prio
     end subroutine _qrm_insert_asm_blk_task_c
  end interface

  interface
     subroutine _qrm_insert_clean_task_c(qrm_dscr, qrm_mat, fnum, sym_handle, prio) bind(c)
       use iso_c_binding
       use qrm_dscr_mod
       type(qrm_dscr_type)           :: qrm_dscr
       type(c_ptr), value            :: qrm_mat, sym_handle
       integer(c_int), value         :: fnum, prio
     end subroutine _qrm_insert_clean_task_c
  end interface
    
  interface
     subroutine _qrm_insert_apply_task_c(qrm_dscr, qrm_mat, b, &
          & front_handle,                            &
          & node_handle,                             &
          & kids_handles,                            &
          & work_handle,                             &
          & transp, inode, nc, prio) bind(c)
       use iso_c_binding
       use qrm_dscr_mod
       type(qrm_dscr_type)           :: qrm_dscr
       type(c_ptr), value            :: qrm_mat, b, node_handle
       type(c_ptr), value            :: work_handle, front_handle
       type(c_ptr)                   :: kids_handles(*)
       character(kind=c_char), value :: transp
       integer(c_int), value         :: nc, inode, prio
     end subroutine _qrm_insert_apply_task_c
  end interface

  
  interface
     subroutine _qrm_insert_solve_task_c(qrm_dscr, qrm_mat, &
          & b, x,                                           &
          & front_handle,                                   &
          & x_node_handle,                                  &
          & b_node_handle,                                  &
          & kids_handles,                                   &
          & transp, inode, nc, prio) bind(c)
       use iso_c_binding
       use qrm_dscr_mod
       type(qrm_dscr_type)           :: qrm_dscr
       type(c_ptr), value            :: qrm_mat, b, x, x_node_handle
       type(c_ptr), value            :: b_node_handle, front_handle
       type(c_ptr)                   :: kids_handles(*)
       character(kind=c_char), value :: transp
       integer(c_int), value         :: nc, inode, prio
     end subroutine _qrm_insert_solve_task_c
  end interface

  interface
     subroutine _qrm_insert_analyse_task_c(qrm_dscr, qrm_mat, &
          & adata_handle,                           &
          & transp) bind(c)
       use qrm_dscr_mod
       use iso_c_binding
       type(qrm_dscr_type)           :: qrm_dscr
       type(c_ptr), value            :: qrm_mat, adata_handle
       character(kind=c_char), value :: transp
     end subroutine _qrm_insert_analyse_task_c
  end interface

  
  interface
     subroutine _qrm_starpu_unpack_args_geqrt(cl_arg, bk, bi, mbs, nbs, ib, cstair) bind(C)
       use iso_c_binding
       type(c_ptr), value            :: cl_arg, bk, bi, mbs, nbs, ib
       type(c_ptr), value            :: cstair
     end subroutine _qrm_starpu_unpack_args_geqrt
  end interface

  interface
     subroutine _qrm_starpu_unpack_args_gemqrt(cl_arg, bk, bi, bj, mbs, nbs, ib, cstair) bind(C)
       use iso_c_binding
       type(c_ptr), value            :: cl_arg, bk, bi, bj, mbs, nbs, ib
       type(c_ptr), value            :: cstair
     end subroutine _qrm_starpu_unpack_args_gemqrt
  end interface

  interface
     subroutine _qrm_starpu_unpack_args_tpqrt(cl_arg, bk, bi, bl, mbs, nbs, ib, ts, cstair) bind(C)
       use iso_c_binding
       type(c_ptr), value            :: cl_arg, bk, bi, bl, mbs, nbs, ib, ts
       type(c_ptr), value            :: cstair
     end subroutine _qrm_starpu_unpack_args_tpqrt
  end interface

  interface
     subroutine _qrm_starpu_unpack_args_tpmqrt(cl_arg, bk, bi, bl, bj, mbs, nbs, ib, ts, cstair) bind(C)
       use iso_c_binding
       type(c_ptr), value            :: cl_arg, bk, bi, bl, bj, mbs, nbs, ib, ts
       type(c_ptr), value            :: cstair
     end subroutine _qrm_starpu_unpack_args_tpmqrt
  end interface

  interface
     subroutine _qrm_starpu_unpack_args_subtree(cl_arg, qrm_dscr_c, qrm_mat_c, rootnum) bind(C)
       use iso_c_binding
       type(c_ptr), value            :: cl_arg, rootnum, qrm_mat_c, qrm_dscr_c
     end subroutine _qrm_starpu_unpack_args_subtree
  end interface

  interface
     subroutine _qrm_starpu_unpack_args_init(cl_arg, qrm_dscr_c, qrm_mat_c, fnum) bind(C)
       use iso_c_binding
       type(c_ptr), value            :: cl_arg, fnum, qrm_mat_c, qrm_dscr_c
     end subroutine _qrm_starpu_unpack_args_init
  end interface

  interface
     subroutine _qrm_starpu_unpack_args_asm(cl_arg, qrm_mat_c, fnum, br, bc) bind(C)
       use iso_c_binding
       type(c_ptr), value            :: cl_arg, fnum, qrm_mat_c, br, bc
     end subroutine _qrm_starpu_unpack_args_asm
  end interface

  interface
     subroutine _qrm_starpu_unpack_args_clean(cl_arg, qrm_dscr, qrm_mat_c, fnum) bind(C)
       use iso_c_binding
       type(c_ptr), value            :: cl_arg, fnum, qrm_mat_c, qrm_dscr
     end subroutine _qrm_starpu_unpack_args_clean
  end interface

  interface
     subroutine _qrm_starpu_unpack_args_apply(cl_arg, qrm_dscr_c, qrm_mat_c, b, inode, transp) bind(C)
       use iso_c_binding
       type(c_ptr), value            :: cl_arg, qrm_mat_c, b, inode, transp, qrm_dscr_c
     end subroutine _qrm_starpu_unpack_args_apply
  end interface

  interface
     subroutine _qrm_starpu_unpack_args_solve(cl_arg, qrm_dscr_c, qrm_mat_c, b, x, inode, transp) bind(C)
       use iso_c_binding
       type(c_ptr), value            :: cl_arg, qrm_dscr_c, qrm_mat_c, b, x, inode, transp
     end subroutine _qrm_starpu_unpack_args_solve
  end interface

  interface
     subroutine _qrm_starpu_unpack_args_analyse(cl_arg, qrm_dscr_c, qrm_mat_c, transp) bind(C)
       use iso_c_binding
       type(c_ptr), value            :: cl_arg, qrm_dscr_c, qrm_mat_c, transp
     end subroutine _qrm_starpu_unpack_args_analyse
  end interface


end module _qrm_starpu_codelets_mod
