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

module qrm_dscr_mod
  use iso_c_binding
  use qrm_error_mod
  use qrm_pthread_mod
  
  type, bind(c) :: qrm_dscr_type
     !> @brief The error status
     integer(c_int)          :: err_status
#if defined (have_starpu)
     !> @brief A lock to prevent simultaneous writing of the error
     type(qrm_pthread_mutex) :: mutex
     integer(c_int)          :: ctx
#endif
  end type qrm_dscr_type

  interface qrm_barrier
     module procedure qrm_barrier, qrm_barrier_dscr
  end interface qrm_barrier

  integer :: qrm_nthreads
  
contains

  subroutine qrm_dscr_init(qrm_dscr)
#if defined (have_starpu)
    use starpu_f_mod
#endif
    use qrm_string_mod
    implicit none

    type(qrm_dscr_type)  :: qrm_dscr
#if defined (have_starpu)
    integer              :: nworkers, i
    integer, allocatable :: workers(:)
#endif

    qrm_dscr%err_status = 0

#if defined (have_starpu)
    call qrm_pthread_mutex_init(qrm_dscr%mutex)

    nworkers = starpu_f_worker_get_count()
    allocate(workers(nworkers))
    do i=1, nworkers
       workers(i) = i-1
    end do
    qrm_dscr%ctx = starpu_f_sched_ctx_create(workers, nworkers, qrm_f2c_str("ctx1"))
    deallocate(workers)
    
#endif

    return

  end subroutine qrm_dscr_init

  subroutine qrm_dscr_destroy(qrm_dscr)

#if defined (have_starpu)
    use starpu_f_mod
#endif
    implicit none
    
    type(qrm_dscr_type) :: qrm_dscr

#if defined (have_starpu)
    call qrm_pthread_mutex_destroy(qrm_dscr%mutex)
    call starpu_f_sched_ctx_delete(qrm_dscr%ctx)
#endif
    
    qrm_dscr%err_status = 0
    
    return

  end subroutine qrm_dscr_destroy

  subroutine qrm_init(nthreads, info)
    use qrm_memhandling_mod
#if defined (have_starpu)
    use starpu_f_mod
#endif
    implicit none

    integer, optional           :: nthreads, info

    integer                     :: ret, ierr
    character(LEN=10)           :: str

    ! error management
    integer                     :: err
    character(len=*), parameter :: name='qrm_init'
    
    err = 0

    if(present(nthreads)) then
       qrm_nthreads = nthreads
    else
       call get_environment_variable(name="QRM_NUM_THREADS",value=str, status=ierr)
       if(ierr .eq. 1) then
          qrm_nthreads = 1
       else
          read(str,*)qrm_nthreads
       end if
    end if
        
#if defined (have_starpu)
    ret = starpu_f_init(max(qrm_nthreads,1))
    if(ret.ne.0) then
       err = 31
       call qrm_error_print(err, name)
    end if
#endif

    qrm_tot_mem = 0
    qrm_max_mem = 0
    
    if(present(info)) info = err
    return
  end subroutine qrm_init


  subroutine qrm_finalize()
    use qrm_memhandling_mod
#if defined (have_starpu)
    use starpu_f_mod
#endif

#if defined (have_starpu)
    call starpu_f_shutdown()
#endif
    
    return
  end subroutine qrm_finalize

  subroutine qrm_barrier()
#if defined (have_starpu)
    use starpu_f_mod
#endif
    implicit none    

#if defined (have_starpu)
    call starpu_f_task_wait_for_all()
#endif
    
    return
  end subroutine qrm_barrier

  subroutine qrm_barrier_dscr(qrm_dscr)
#if defined (have_starpu)
    use starpu_f_mod 
#endif
    implicit none    

    type(qrm_dscr_type) :: qrm_dscr

    
#if defined (have_starpu)
    call starpu_f_task_wait_for_all_in_ctx(qrm_dscr%ctx)
#endif
    
    return
  end subroutine qrm_barrier_dscr


  subroutine qrm_status_set(qrm_dscr, code, where, ied, aed)
    use qrm_error_mod
    use qrm_const_mod
    implicit none
    
    type(qrm_dscr_type)             :: qrm_dscr
    integer                         :: code
    character(len=*)                :: where
    integer, optional               :: ied(:)
    character(len=*), optional      :: aed

#if defined (have_starpu)
    call qrm_pthread_mutex_lock(qrm_dscr%mutex)
#endif
    if(qrm_dscr%err_status.eq.0) then
       qrm_dscr%err_status = code
       call qrm_error_print(code, where, ied, aed)
    end if
#if defined (have_starpu)
    call qrm_pthread_mutex_unlock(qrm_dscr%mutex)
#endif
    return
  end subroutine qrm_status_set
  
  function qrm_get_num_threads(qrm_dscr) result(nworkers)
    use qrm_error_mod
    use qrm_const_mod
#if defined (have_starpu)
    use starpu_f_mod
#endif
    implicit none
    
    type(qrm_dscr_type), optional :: qrm_dscr
    integer                       :: nworkers
    
    if(present(qrm_dscr)) then
#if defined (have_starpu)
       nworkers = starpu_f_sched_ctx_get_nworkers(qrm_dscr%ctx)
#else             
       nworkers = 1
#endif
    else
#if defined (have_starpu)
       nworkers = starpu_f_cpu_worker_get_count()
#else
       nworkers = 1
#endif
    end if
    
    return
  end function qrm_get_num_threads
  
  
end module qrm_dscr_mod
