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

module qrm_memhandling_mod
  use iso_c_binding
  use qrm_error_mod
  use qrm_const_mod
  use qrm_pthread_mod

  integer, parameter :: qrm_mem_kind=c_int64_t

  !>  a counter to keep track of the currently allocated memory, per thread
  integer(kind=qrm_mem_kind), save :: qrm_tot_mem=0

  !>  a counter to keep track of the peak memory, per thread
  integer(kind=qrm_mem_kind), save :: qrm_max_mem=0

  integer(kind=qrm_mem_kind), parameter :: qrm_sizeof_i_  = 4
  integer(kind=qrm_mem_kind), parameter :: qrm_sizeof_i8_ = 4
  integer(kind=qrm_mem_kind), parameter :: qrm_sizeof_s_  = 4
  integer(kind=qrm_mem_kind), parameter :: qrm_sizeof_d_  = 8
  integer(kind=qrm_mem_kind), parameter :: qrm_sizeof_c_  = 8
  integer(kind=qrm_mem_kind), parameter :: qrm_sizeof_z_  = 16


  type qrm_ma_type
     integer(kind=qrm_mem_kind) :: avail=0, peak, cons
     ! integer(kind=qrm_mem_kind) :: peak, cons
     type(qrm_pthread_mutex)    :: mutex
     type(qrm_pthread_cond)     :: cond
  end type qrm_ma_type
  
contains

  !> @brief This routine has to be called at the beginning of a
  !! factorization and initializes whatever is needed for the
  !! memory-aware tasks submission
  subroutine qrm_facto_mem_init(ma, available_mem)
    use qrm_pthread_mod
    implicit none
    type(qrm_ma_type)          :: ma
    integer(kind=qrm_mem_kind) :: available_mem
    integer :: info

    ma%avail = available_mem
    ! ma%peak  = 0
    ! ma%cons  = 0
    
    call qrm_pthread_mutex_init(ma%mutex, info)
    call qrm_pthread_cond_init(ma%cond, info)
    
    return
  end subroutine qrm_facto_mem_init
    
  subroutine qrm_facto_mem_finalize(ma)
    use qrm_pthread_mod
    implicit none
    type(qrm_ma_type)  :: ma
    integer            :: info

    call qrm_pthread_mutex_destroy(ma%mutex, info)
    call qrm_pthread_cond_destroy(ma%cond, info)

    return
  end subroutine qrm_facto_mem_finalize

  
  subroutine qrm_facto_mem_get(ma, mem)
    use qrm_pthread_mod
    implicit none
    type(qrm_ma_type)          :: ma
    integer(kind=qrm_mem_kind) :: mem

    call qrm_pthread_mutex_lock(ma%mutex)
    do while(ma%avail .lt. mem)
       call qrm_pthread_cond_wait(ma%cond, &
            & ma%mutex)
    end do

    ma%avail = ma%avail - mem
    ! ma%cons  = ma%cons + mem
    ! ma%peak  = max(ma%peak, ma%cons)

    call qrm_pthread_cond_signal(ma%cond)
    
    call qrm_pthread_mutex_unlock(ma%mutex)
    
    return
  end subroutine qrm_facto_mem_get
  
  
  !> @brief updates memory statistics
  !!
  !! @param[in] n the amount of memory to be added (can be negative
  !!            for deallocations)
  subroutine qrm_mem_upd(n)
    use qrm_pthread_mod
    implicit none
    
    integer(kind=qrm_mem_kind) :: n

    call qrm_atomic_add(qrm_tot_mem, n)
    call qrm_atomic_max(qrm_max_mem, qrm_tot_mem)

    return
  end subroutine qrm_mem_upd
  
end module qrm_memhandling_mod
