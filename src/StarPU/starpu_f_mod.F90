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


module starpu_f_mod
  use iso_c_binding

  interface starpu_f_shutdown
     subroutine starpu_shutdown() bind(c)
     end subroutine starpu_shutdown
  end interface starpu_f_shutdown

  interface starpu_f_pause
     subroutine starpu_pause() bind(c)
     end subroutine starpu_pause
  end interface starpu_f_pause

  interface starpu_f_resume
     subroutine starpu_resume() bind(c)
     end subroutine starpu_resume
  end interface starpu_f_resume

  interface starpu_f_init
     function starpu_f_init_c(ncpus) bind(c)
       use iso_c_binding
       integer(c_int)        :: starpu_f_init_c
       integer(c_int), value :: ncpus 
     end function starpu_f_init_c
     function starpu_init(conf) bind(c)
       use iso_c_binding
       integer(c_int)     :: starpu_init
       type(c_ptr), value :: conf
     end function starpu_init
  end interface starpu_f_init

  ! contexts
  interface starpu_f_sched_ctx_create
     function starpu_f_sched_ctx_create_c(workers, nworkers, ctxname) bind(c)
       use iso_c_binding
       integer(c_int)         :: starpu_f_sched_ctx_create_c
       integer(c_int)         :: workers(*)
       integer(c_int), value  :: nworkers
       character(kind=c_char) :: ctxname(*)
     end function starpu_f_sched_ctx_create_c
  end interface starpu_f_sched_ctx_create

  interface starpu_f_sched_ctx_set_context
     subroutine starpu_f_sched_ctx_set_context_c(ctx) bind(c)
       use iso_c_binding
       integer(c_int)         :: ctx
     end subroutine starpu_f_sched_ctx_set_context_c
  end interface starpu_f_sched_ctx_set_context
  
  interface starpu_f_sched_ctx_display_workers
     subroutine starpu_f_sched_ctx_display_workers_c(ctx) bind(c)
       use iso_c_binding
       integer(c_int), value      :: ctx
     end subroutine starpu_f_sched_ctx_display_workers_c
  end interface starpu_f_sched_ctx_display_workers
  
  interface starpu_f_sched_ctx_delete
     subroutine starpu_f_sched_ctx_delete_c(ctx) bind(c)
       use iso_c_binding
       integer(c_int), value      :: ctx
     end subroutine starpu_f_sched_ctx_delete_c
  end interface starpu_f_sched_ctx_delete
  
  interface starpu_f_task_wait_for_all_in_ctx
     subroutine starpu_f_task_wait_for_all_in_ctx_c(ctx) bind(c)
       use iso_c_binding
       integer(c_int), value      :: ctx
     end subroutine starpu_f_task_wait_for_all_in_ctx_c
  end interface starpu_f_task_wait_for_all_in_ctx
  
  ! workers

  ! return the number of workers
  interface starpu_f_cpu_worker_get_count
     function starpu_cpu_worker_get_count() bind(c)
       use iso_c_binding
       integer(c_int) :: starpu_cpu_worker_get_count
     end function starpu_cpu_worker_get_count
  end interface starpu_f_cpu_worker_get_count

  interface starpu_f_sched_ctx_get_nworkers
     function starpu_sched_ctx_get_nworkers (ctx) bind(c)
       use iso_c_binding
       integer(c_int)        :: starpu_sched_ctx_get_nworkers
       integer(c_int), value :: ctx
     end function starpu_sched_ctx_get_nworkers
  end interface starpu_f_sched_ctx_get_nworkers
  
  ! return the number of workers
  interface starpu_f_worker_get_count
     function starpu_worker_get_count() bind(c)
       use iso_c_binding
       integer(c_int) :: starpu_worker_get_count
     end function starpu_worker_get_count
  end interface starpu_f_worker_get_count
  
  ! return the id of the current worker
  interface starpu_f_worker_get_id
     function starpu_worker_get_id() bind(c)
       use iso_c_binding
       integer(c_int) :: starpu_worker_get_id
     end function starpu_worker_get_id
  end interface starpu_f_worker_get_id

  ! data interfaces 

  interface starpu_f_void_data_register
     subroutine starpu_void_data_register(handle) bind(c)
       use iso_c_binding
       type(c_ptr)    :: handle
     end subroutine starpu_void_data_register
  end interface starpu_f_void_data_register

  interface starpu_f_data_unregister
     subroutine starpu_data_unregister(handle) bind(c)
       use iso_c_binding
       type(c_ptr), value            :: handle
     end subroutine starpu_data_unregister
  end interface starpu_f_data_unregister
  
  interface starpu_f_data_unregister_no_coherency
     subroutine starpu_data_unregister_no_coherency(handle) bind(c)
       use iso_c_binding
       type(c_ptr), value            :: handle
     end subroutine starpu_data_unregister_no_coherency
  end interface starpu_f_data_unregister_no_coherency
  

  interface starpu_f_data_acquire_read
     subroutine starpu_f_data_acquire_read(handle) bind(c)
       use iso_c_binding
       type(c_ptr), value            :: handle
     end subroutine starpu_f_data_acquire_read
  end interface starpu_f_data_acquire_read
  
  interface starpu_f_data_release
     subroutine starpu_data_release(handle) bind(c)
       use iso_c_binding
       type(c_ptr), value            :: handle
     end subroutine starpu_data_release
  end interface starpu_f_data_release
  
  
  ! tasks

  interface starpu_f_task_wait_for_all
     subroutine starpu_task_wait_for_all() bind(c)
     end subroutine starpu_task_wait_for_all
  end interface starpu_f_task_wait_for_all

  interface starpu_f_get_buffer
     subroutine starpu_f_get_buffer(buffers, num, a, m, n, lda) bind(C)
       use iso_c_binding
       type(c_ptr), value    :: a
       integer(c_int), value :: num
       type(c_ptr), value    :: m, n, lda, buffers
     end subroutine starpu_f_get_buffer
  end interface starpu_f_get_buffer



  interface starpu_f_matrix_data_register
     module procedure starpu_f_1dsmatrix_data_register, starpu_f_2dsmatrix_data_register
     module procedure starpu_f_1ddmatrix_data_register, starpu_f_2ddmatrix_data_register
     module procedure starpu_f_1dcmatrix_data_register, starpu_f_2dcmatrix_data_register
     module procedure starpu_f_1dzmatrix_data_register, starpu_f_2dzmatrix_data_register
     subroutine starpu_matrix_data_register(handle, host, a, ld, nx, ny, sizeof) bind(c)
       use iso_c_binding
       type(c_ptr)              :: handle
       type(c_ptr), value       :: a
       integer(c_int), value    :: host, ld, nx, ny
       integer(c_size_t), value :: sizeof
     end subroutine starpu_matrix_data_register
  end interface starpu_f_matrix_data_register


  ! interface
  ! end interface

  interface starpu_f_data_unregister_submit
     subroutine starpu_data_unregister_submit(handle) bind(c)
       use iso_c_binding
       type(c_ptr), value :: handle
     end subroutine starpu_data_unregister_submit
  end interface starpu_f_data_unregister_submit
  

  
  interface starpu_f_fxt_start_profiling
     subroutine starpu_fxt_start_profiling() bind(c)
     end subroutine starpu_fxt_start_profiling
  end interface starpu_f_fxt_start_profiling
  interface starpu_f_fxt_stop_profiling
     subroutine starpu_fxt_stop_profiling() bind(c)
     end subroutine starpu_fxt_stop_profiling
  end interface starpu_f_fxt_stop_profiling

contains

  subroutine starpu_f_1dsmatrix_data_register(handle, host, a, ld, m, n)
    implicit none
    type(c_ptr)       :: handle
    real, target      :: a(:)
    integer           :: host
    integer, optional :: ld, m, n
    integer           :: ild, im, in
    real, pointer     :: ptr
    
    ptr => a(1)
    
    if(present(ld)) then
       ild = ld
    else
       ild = size(a,1)
    end if

    if(present(m)) then
       im = m
    else
       im = size(a,1)
    end if

    if(present(n)) then
       in = n
    else
       in = 1
    end if

    ! 4 is hardcoded. should make something better and more portable 
    call starpu_matrix_data_register(handle, host, c_loc(ptr), ild, im, in, int(4,kind=c_size_t))

    return
  end subroutine starpu_f_1dsmatrix_data_register


  subroutine starpu_f_2dsmatrix_data_register(handle, host, a, ld, m, n)
    implicit none
    type(c_ptr)       :: handle
    real, target      :: a(:,:)
    integer           :: host
    integer, optional :: ld, m, n
    integer           :: ild, im, in
    real, pointer     :: ptr

    ptr => a(1,1)
    
    if(present(ld)) then
       ild = ld
    else
       ild = size(a,1)
    end if

    if(present(m)) then
       im = m
    else
       im = size(a,1)
    end if

    if(present(n)) then
       in = n
    else
       in = size(a,2)
    end if

    ! 4 is hardcoded. should make something better and more portable 
    call starpu_matrix_data_register(handle, host, c_loc(ptr), ild, im, in, int(4,kind=c_size_t))

    return
  end subroutine starpu_f_2dsmatrix_data_register

  subroutine starpu_f_1ddmatrix_data_register(handle, host, a, ld, m, n)
    implicit none
    type(c_ptr)               :: handle
    real(kind(1.d0)), target  :: a(:)
    integer                   :: host
    integer, optional         :: ld, m, n
    integer                   :: ild, im, in
    real(kind(1.d0)), pointer :: ptr

    ptr => a(1)
    
    if(present(ld)) then
       ild = ld
    else
       ild = size(a,1)
    end if

    if(present(m)) then
       im = m
    else
       im = size(a,1)
    end if

    if(present(n)) then
       in = n
    else
       in = 1
    end if

    ! 4 is hardcoded. should make something better and more portable 
    call starpu_matrix_data_register(handle, host, c_loc(ptr), ild, im, in, int(8,kind=c_size_t))

    return
  end subroutine starpu_f_1ddmatrix_data_register


  subroutine starpu_f_2ddmatrix_data_register(handle, host, a, ld, m, n)
    implicit none
    type(c_ptr)               :: handle
    real(kind(1.d0)), target  :: a(:,:)
    integer                   :: host
    integer, optional         :: ld, m, n
    integer                   :: ild, im, in
    real(kind(1.d0)), pointer :: ptr

    ptr => a(1,1)
    
    if(present(ld)) then
       ild = ld
    else
       ild = size(a,1)
    end if

    if(present(m)) then
       im = m
    else
       im = size(a,1)
    end if

    if(present(n)) then
       in = n
    else
       in = size(a,2)
    end if

    ! 4 is hardcoded. should make something better and more portable 
    call starpu_matrix_data_register(handle, host, c_loc(ptr), ild, im, in, int(8,kind=c_size_t))

    return
  end subroutine starpu_f_2ddmatrix_data_register

  subroutine starpu_f_1dcmatrix_data_register(handle, host, a, ld, m, n)
    implicit none
    type(c_ptr)       :: handle
    complex, target   :: a(:)
    integer           :: host
    integer, optional :: ld, m, n
    integer           :: ild, im, in
    complex, pointer  :: ptr

    ptr => a(1)
    
    if(present(ld)) then
       ild = ld
    else
       ild = size(a,1)
    end if

    if(present(m)) then
       im = m
    else
       im = size(a,1)
    end if

    if(present(n)) then
       in = n
    else
       in = 1
    end if

    ! 4 is hardcoded. should make something better and more portable 
    call starpu_matrix_data_register(handle, host, c_loc(ptr), ild, im, in, int(8,kind=c_size_t))

    return
  end subroutine starpu_f_1dcmatrix_data_register



  subroutine starpu_f_2dcmatrix_data_register(handle, host, a, ld, m, n)
    implicit none
    type(c_ptr)       :: handle
    complex, target   :: a(:,:)
    integer           :: host
    integer, optional :: ld, m, n
    integer           :: ild, im, in
    complex, pointer  :: ptr

    ptr => a(1,1)
    
    if(present(ld)) then
       ild = ld
    else
       ild = size(a,1)
    end if

    if(present(m)) then
       im = m
    else
       im = size(a,1)
    end if

    if(present(n)) then
       in = n
    else
       in = size(a,2)
    end if

    ! 4 is hardcoded. should make something better and more portable 
    call starpu_matrix_data_register(handle, host, c_loc(ptr), ild, im, in, int(8,kind=c_size_t))

    return
  end subroutine starpu_f_2dcmatrix_data_register

  subroutine starpu_f_1dzmatrix_data_register(handle, host, a, ld, m, n)
    implicit none
    type(c_ptr)                  :: handle
    complex(kind(1.d0)), target  :: a(:)
    integer                      :: host
    integer, optional            :: ld, m, n
    integer                      :: ild, im, in
    complex(kind(1.d0)), pointer :: ptr

    ptr => a(1)
    
    if(present(ld)) then
       ild = ld
    else
       ild = size(a,1)
    end if

    if(present(m)) then
       im = m
    else
       im = size(a,1)
    end if

    if(present(n)) then
       in = n
    else
       in = 1
    end if

    ! 4 is hardcoded. should make something better and more portable 
    call starpu_matrix_data_register(handle, host, c_loc(ptr), ild, im, in, int(16,kind=c_size_t))

    return
  end subroutine starpu_f_1dzmatrix_data_register


  subroutine starpu_f_2dzmatrix_data_register(handle, host, a, ld, m, n)
    implicit none
    type(c_ptr)                  :: handle
    complex(kind(1.d0)), target  :: a(:,:)
    integer                      :: host
    integer, optional            :: ld, m, n
    integer                      :: ild, im, in
    complex(kind(1.d0)), pointer :: ptr

    ptr => a(1,1)
    
    if(present(ld)) then
       ild = ld
    else
       ild = size(a,1)
    end if

    if(present(m)) then
       im = m
    else
       im = size(a,1)
    end if

    if(present(n)) then
       in = n
    else
       in = size(a,2)
    end if

    ! 4 is hardcoded. should make something better and more portable 
    call starpu_matrix_data_register(handle, host, c_loc(ptr), ild, im, in, int(16,kind=c_size_t))

    return
  end subroutine starpu_f_2dzmatrix_data_register

end module starpu_f_mod
