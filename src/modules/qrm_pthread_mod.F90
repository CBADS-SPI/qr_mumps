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


module qrm_pthread_mod
  use iso_c_binding

  ! private

  ! public :: qrm_atomic_max
  ! public :: qrm_atomic_add
  ! public :: qrm_pthread_mutex_init
  ! public :: qrm_pthread_mutex_destroy
  ! public :: qrm_pthread_mutex_lock
  ! public :: qrm_pthread_mutex_unlock
  ! public :: qrm_pthread_cond_init
  ! public :: qrm_pthread_cond_destroy
  ! public :: qrm_pthread_cond_wait
  ! public :: qrm_pthread_cond_signal
  
  ! not much to do with pthreads but still useful for threading

  ! atomic max of  variable and a value
  interface qrm_atomic_max
     subroutine qrm_atomic_max_int64_t(ptr, val) bind(C)
       use iso_c_binding
       integer(kind=c_int64_t)        :: ptr
       integer(kind=c_int64_t), value :: val
     end subroutine qrm_atomic_max_int64_t
  end interface qrm_atomic_max
  
  ! atomic increment of a variable
  interface qrm_atomic_add
     subroutine qrm_atomic_add_int64_t(ptr, val) bind(C)
       use iso_c_binding
       integer(kind=c_int64_t)        :: ptr
       integer(kind=c_int64_t), value :: val
     end subroutine qrm_atomic_add_int64_t
  end interface qrm_atomic_add

  ! type mutex 
  type, bind(c) :: qrm_pthread_mutex
     type(c_ptr)                      :: m
  end type qrm_pthread_mutex

  ! type cond 
  type, public :: qrm_pthread_cond
     ! private
     type(c_ptr)                      :: c
  end type qrm_pthread_cond


  ! allocate mutex
  interface qrm_alloc_pthread_mutex
     subroutine qrm_alloc_pthread_mutex_c(mutex) bind(c)
       use iso_c_binding
       type(c_ptr)                    :: mutex
     end subroutine qrm_alloc_pthread_mutex_c
  end interface qrm_alloc_pthread_mutex

  ! allocate condition
  interface qrm_alloc_pthread_cond
     subroutine qrm_alloc_pthread_cond_c(cond) bind(c)
       use iso_c_binding
       type(c_ptr)                    :: cond
     end subroutine qrm_alloc_pthread_cond_c
  end interface qrm_alloc_pthread_cond

  ! allocate mutex
  interface qrm_dealloc_pthread_mutex
     subroutine qrm_dealloc_pthread_mutex_c(mutex) bind(c)
       use iso_c_binding
       type(c_ptr)                    :: mutex
     end subroutine qrm_dealloc_pthread_mutex_c
  end interface qrm_dealloc_pthread_mutex

  ! allocate condition
  interface qrm_dealloc_pthread_cond
     subroutine qrm_dealloc_pthread_cond_c(cond) bind(c)
       use iso_c_binding
       type(c_ptr)                    :: cond
     end subroutine qrm_dealloc_pthread_cond_c
  end interface qrm_dealloc_pthread_cond


  
contains

  ! initialize mutex
  subroutine qrm_pthread_mutex_init(mutex, info)
    use iso_c_binding
    implicit none
    type(qrm_pthread_mutex)           :: mutex
    integer, optional                 :: info

    interface
       function pthread_mutex_init_c(mutex, attr) bind(c)
         use iso_c_binding
         integer(c_int)               :: pthread_mutex_init_c
         type(c_ptr), value           :: mutex
         type(c_ptr), value           :: attr
       end function pthread_mutex_init_c
    end interface

    integer                           :: err
    
    call qrm_alloc_pthread_mutex(mutex%m)
    err = pthread_mutex_init_c(mutex%m, c_null_ptr)
    if(present(info)) info = err    
    return

  end subroutine qrm_pthread_mutex_init

  ! check if a mutex is initialized
  ! function qrm_pthread_mutex_isinit(mutex)
    ! use iso_c_binding
    ! implicit none
    ! type(qrm_pthread_mutex)           :: mutex
    ! logical                           :: qrm_pthread_mutex_isinit
    ! qrm_pthread_mutex_isinit=c_associated(mutex%m)
    ! return
  ! end function qrm_pthread_mutex_isinit


  ! destroy mutex
  subroutine qrm_pthread_mutex_destroy(mutex, info)
    implicit none

    type(qrm_pthread_mutex)           :: mutex
    integer, optional                 :: info
    
    interface
       function pthread_mutex_destroy(mutex) bind(c)
         use iso_c_binding
         integer(c_int)               :: pthread_mutex_destroy
         type(c_ptr), value           :: mutex
       end function pthread_mutex_destroy
    end interface

    integer                           :: err
    
    err = pthread_mutex_destroy(mutex%m)
    call qrm_dealloc_pthread_mutex(mutex%m)
    if(present(info)) info=err
    return
    
  end subroutine qrm_pthread_mutex_destroy

  
    
  ! initialize condition variable
  subroutine qrm_pthread_cond_init(cond, info)
    use iso_c_binding
    implicit none
    type(qrm_pthread_cond)            :: cond
    integer, optional                 :: info

    interface
       function pthread_cond_init_c(cond, attr) bind(c)
         use iso_c_binding
         integer(c_int)               :: pthread_cond_init_c
         type(c_ptr), value           :: cond
         type(c_ptr), value           :: attr
       end function pthread_cond_init_c
    end interface

    integer                           :: err
    
    call qrm_alloc_pthread_cond(cond%c)
    err = pthread_cond_init_c(cond%c, c_null_ptr)
    if(present(info)) info=err
    return

  end subroutine qrm_pthread_cond_init


  ! check if a condition variable is initialized
  ! function qrm_pthread_cond_isinit(cond)
    ! use iso_c_binding
    ! implicit none
    ! type(qrm_pthread_cond)            :: cond
    ! logical                           :: qrm_pthread_cond_isinit
    ! qrm_pthread_cond_isinit=C_ASSOCIATED(cond%c)
    ! return
  ! end function qrm_pthread_cond_isinit


  ! destroy mutex
  subroutine qrm_pthread_cond_destroy(cond, info)
    implicit none
    
    type(qrm_pthread_cond)            :: cond
    integer, optional                 :: info
    
    ! destroy condition variable
    interface 
       function pthread_cond_destroy(cond) bind(c)
         use iso_c_binding
         integer(c_int)               :: pthread_cond_destroy
         type(c_ptr), value           :: cond
       end function pthread_cond_destroy
    end interface

    integer                           :: err
    
    err = pthread_cond_destroy(cond%c)
    call qrm_dealloc_pthread_cond(cond%c)
    if(present(info)) info=err
    return
    
  end subroutine qrm_pthread_cond_destroy

  

  
  ! wrapper for pthread_mutex_lock function 
  subroutine qrm_pthread_mutex_lock(mutex, info)
    implicit none

    type(qrm_pthread_mutex)           :: mutex
    integer, optional                 :: info
    
    interface 
       function pthread_mutex_lock(mutex) bind(c)
         use iso_c_binding
         integer(c_int)               :: pthread_mutex_lock
         type(c_ptr), value           :: mutex
       end function pthread_mutex_lock
    end interface

    integer                           :: err
    
    err = pthread_mutex_lock(mutex%m)
    if(present(info)) info=err
    return
    
  end subroutine qrm_pthread_mutex_lock

  ! wrapper for pthread_mutex_unlock function 
  subroutine qrm_pthread_mutex_unlock(mutex, info)
    implicit none
    
    type(qrm_pthread_mutex)           :: mutex
    integer, optional                 :: info
        
    ! unlock mutex
    interface 
       function pthread_mutex_unlock(mutex) bind(c)
         use iso_c_binding
         integer(c_int)               :: pthread_mutex_unlock
         type(c_ptr), value           :: mutex
       end function pthread_mutex_unlock
    end interface
    
    integer                           :: err
    
    err = pthread_mutex_unlock(mutex%m)
    if(present(info)) info=err
    return

  end subroutine qrm_pthread_mutex_unlock

  ! wrapper for pthread_cond_wait function
  subroutine qrm_pthread_cond_wait(cond, mutex, info)
    implicit none
    type(qrm_pthread_cond)            :: cond
    type(qrm_pthread_mutex)           :: mutex
    integer, optional                 :: info

    interface
       function pthread_cond_wait_c(cond, mutex) bind(c)
         use iso_c_binding
         integer(c_int)               :: pthread_cond_wait_c
         type(c_ptr), value           :: cond
         type(c_ptr), value           :: mutex       
       end function pthread_cond_wait_c    
    end interface

    integer                           :: err
    
    err = pthread_cond_wait_c(cond%c, mutex%m)
    if(present(info)) info=err
    return 

  end subroutine qrm_pthread_cond_wait

  ! wrapper for pthread_cond_wait function
  subroutine qrm_pthread_cond_signal(cond, info)
    implicit none
    type(qrm_pthread_cond)            :: cond
    integer, optional                 :: info

    interface
       function pthread_cond_signal_c(cond) bind(c)
         use iso_c_binding
         integer(c_int)               :: pthread_cond_signal_c
         type(c_ptr), value           :: cond    
       end function pthread_cond_signal_c
    end interface

    integer                           :: err
    
    err = pthread_cond_signal_c(cond%c)
    if(present(info)) info=err
    return 

  end subroutine qrm_pthread_cond_signal

  
end module qrm_pthread_mod
