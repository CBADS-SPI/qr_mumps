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
!> @file qrm_test_c.F90
!! This file contains tests for the C interface
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################

function _qrm_test_c(c, m) result(ok)
  use _qrm_testing_mod, protect => _qrm_test_c
  use _qrm_c_interface
  implicit none
  logical :: ok

  integer :: c, m

  integer, parameter :: ncases=15
  logical :: cases(ncases)

  if(c .eq. -1) then
     cases = .true.
  else if (c .le. ncases) then
     cases = .false.
     cases(c) = .true.
  end if

  ok = .true.
  
  if(cases(1 )) ok = case1 (m) .and. ok 
  if(cases(2 )) ok = case2 (m) .and. ok 


  return


contains

  function case1(m) result(ok)
    use _qrm_mod
    use _qrm_c_interface
    implicit none

    logical :: ok
    integer :: m

    interface
       function _qrm_test_solve_c(qrm_spmat_c, b_c, x_c, r_c, eps, err) result(info) bind(c)
         use iso_c_binding
         use _qrm_c_interface
         _qrm_real_fc            :: err
         _qrm_real_fc, value     :: eps
         type(_qrm_spmat_type_c) :: qrm_spmat_c
         type(c_ptr), value      :: b_c, x_c, r_c
         integer(c_int)          :: info
       end function _qrm_test_solve_c
    end interface
    

    type(_qrm_spmat_type), pointer :: qrm_spmat
    integer                        :: info
    _qrm_data, pointer             :: x(:), b(:), r(:)
    _qrm_real                      :: err
    integer                        :: i
    type(_qrm_spmat_type_c)        :: qrm_spmat_c
    type(c_ptr)                    :: b_c, x_c, r_c

    nullify(b) 
    nullify(x) 
    nullify(r) 
    
    if((m.ne.-1) .and. ((m.lt.11).and.(m.gt.size(matrices)))) then
       write(*,'("Matrix ",i2," is not available for this test")')m
       return
    end if

    ok = .true.
    
    ! loop over all the file-matrices
    do i=11, size(matrices)
       if((m.ne.-1) .and. (i.ne.m)) cycle
       info = 0
       
       ! get the data
       qrm_spmat => _qrm_get_test_mat(i)
       if(info.eq.0) call qrm_alloc(b, qrm_spmat%m, info)
       if(info.eq.0) call qrm_alloc(r, qrm_spmat%m, info)
       if(info.eq.0) call qrm_alloc(x, qrm_spmat%n, info)
       if(info.ne.0) goto 9999
       call _xlarnv(2, iseed, size(b), b(1))
       r = b
          
       info = _qrm_spmat_init_c(qrm_spmat_c)
       if(info.ne.0) goto 9999

       qrm_spmat_c%m   = qrm_spmat%m
       qrm_spmat_c%n   = qrm_spmat%n
       qrm_spmat_c%nz  = qrm_spmat%nz

       qrm_spmat_c%irn = c_loc(qrm_spmat%irn(1))
       qrm_spmat_c%jcn = c_loc(qrm_spmat%jcn(1))
       qrm_spmat_c%val = c_loc(qrm_spmat%val(1))
       b_c             = c_loc(b(1))
       x_c             = c_loc(x(1))
       r_c             = c_loc(r(1))
          
       ! solve and get back the error
       info = _qrm_test_solve_c(qrm_spmat_c, b_c, x_c, r_c, eps, err)
       
       ! print message
       call _qrm_prnt_testmesg(3, "c_int", 1, 1, i, (info .eq. 0) .and. (err .lt. eps))
       ok = ok .and. (info .eq. 0) .and. (err .lt. eps)
       
       info = _qrm_spmat_destroy_c(qrm_spmat_c)
       ! put cntl back in its initial state
       call _qrm_spmat_cleanup(qrm_spmat)
       call _qrm_cntl_init(qrm_spmat)

9999   continue
       call qrm_dealloc(b)
       call qrm_dealloc(r)
       call qrm_dealloc(x)

    end do

    return

  end function case1

  
  !> @brief Test C problem solve
  function case2(m) result(ok)
    use _qrm_mod
    use _qrm_c_interface
    implicit none
    logical :: ok
    integer :: m

    interface
       function _qrm_test_methods_c(qrm_spmat_c, b_c, x_c, r_c, eps, err) result(info) bind(c)
         use iso_c_binding
         use _qrm_c_interface
         _qrm_real_fc            :: err
         _qrm_real_fc, value     :: eps
         type(_qrm_spmat_type_c) :: qrm_spmat_c
         type(c_ptr), value      :: b_c, x_c, r_c
         integer(c_int)          :: info
       end function _qrm_test_methods_c
    end interface
    
    type(_qrm_spmat_type), pointer :: qrm_spmat
    integer                        :: info
    _qrm_data, pointer             :: x(:), b(:), r(:)
    _qrm_real                      :: err
    integer                        :: i
    type(_qrm_spmat_type_c)        :: qrm_spmat_c
    type(c_ptr)                    :: b_c, x_c, r_c

    nullify(b) 
    nullify(x) 
    nullify(r) 
    
    if((m.ne.-1) .and. ((m.lt.11).and.(m.gt.size(matrices)))) then
       write(*,'("Matrix ",i2," is not available for this test")')m
       return
    end if

    ok = .true.
    
    ! loop over all the file-matrices
    do i=11, size(matrices)
       if((m.ne.-1) .and. (i.ne.m)) cycle
       info = 0
       
       ! get the data
       qrm_spmat => _qrm_get_test_mat(i)
       if(info.eq.0) call qrm_alloc(b, qrm_spmat%m, info)
       if(info.eq.0) call qrm_alloc(r, qrm_spmat%m, info)
       if(info.eq.0) call qrm_alloc(x, qrm_spmat%n, info)
       if(info.ne.0) goto 9999
       call _xlarnv(2, iseed, size(b), b(1))
       r = b
          
       info = _qrm_spmat_init_c(qrm_spmat_c)
       if(info.ne.0) goto 9999

       qrm_spmat_c%m   = qrm_spmat%m
       qrm_spmat_c%n   = qrm_spmat%n
       qrm_spmat_c%nz  = qrm_spmat%nz

       qrm_spmat_c%irn = c_loc(qrm_spmat%irn(1))
       qrm_spmat_c%jcn = c_loc(qrm_spmat%jcn(1))
       qrm_spmat_c%val = c_loc(qrm_spmat%val(1))
       b_c             = c_loc(b(1))
       x_c             = c_loc(x(1))
       r_c             = c_loc(r(1))
          
       ! solve and get back the error
       info = _qrm_test_methods_c(qrm_spmat_c, b_c, x_c, r_c, eps, err)
       
       ! print message
       call _qrm_prnt_testmesg(3, "c_int", 2, 1, i, (info .eq. 0) .and. (err .lt. eps))
       ok = ok .and. (info .eq. 0) .and. (err .lt. eps)          
       info = _qrm_spmat_destroy_c(qrm_spmat_c)
       ! put cntl back in its initial state
       call _qrm_spmat_cleanup(qrm_spmat)
       call _qrm_cntl_init(qrm_spmat)

9999   continue
       call qrm_dealloc(b)
       call qrm_dealloc(r)
       call qrm_dealloc(x)

    end do

    return

  end function case2
  



end function _qrm_test_c
