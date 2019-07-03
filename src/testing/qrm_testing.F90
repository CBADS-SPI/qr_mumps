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
!> @file qrm_testing.F90
!! This file contains a testing test program to check as many features and lines of
!! code as possible.
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


program _qrm_testing
  use _qrm_mod
  use _qrm_testing_mod
  implicit none

  integer :: t, c, m
  integer :: nargs, i
  integer, parameter :: ntests=7
  logical :: tests(ntests), ok


  call get_args(t,c,m)
  call read_matfile()

  call qrm_set('qrm_dunit',-1)
  call qrm_set('qrm_eunit',-1)
  call qrm_set('qrm_ounit',-1)

  write(*,'("=============================================================================")')
  write(*,'("      _qrm_testing ")')
  write(*,'(" ")')
  write(*,'("Matrices used for the test")')
  do i=11, size(matrices)
     write(*,'(i2," -- ",a30)')i,matrices(i)%mfile
  end do

  write(*,'(" ")')
  if(t .eq. -1) then
     tests = .true.
  else
     tests = .false.
     tests(t) = .true.
  end if

  ok = .true.
  
  call qrm_init()
  if(tests(1)) ok = _qrm_test_err(c)          .and. ok 
  if(tests(2)) ok = _qrm_test_ord(c, m)       .and. ok 
  if(tests(3)) ok = _qrm_test_c(c, m)         .and. ok 
  if(tests(4)) ok = _qrm_test_multiples(c, m) .and. ok 
  if(tests(5)) ok = _qrm_test_facto(c, m)     .and. ok 
  if(tests(6)) ok = _qrm_test_methods(c, m)   .and. ok 
  if(tests(7)) ok = _qrm_test_pipe(c, m)      .and. ok 

  call qrm_finalize()

  write(*,'(" ")')
  if(ok) then
     write(*,'("All tests suceeded")')
     write(*,'("=============================================================================")')
  else
     write(*,'("Some tests failed")')
     write(*,'("=============================================================================")')
     error stop
  end if

  stop
  
contains
  
  subroutine get_args(t, c, m)
    
    integer :: t, c, m
    
    character(len=50) :: str
    integer :: idx, len, i

    t = -1
    m = -1
    c = -1
    
    nargs = command_argument_count()
    if(nargs .gt. 0) then
       call get_command_argument(1,value=str,length=len)
       idx = index(str,'-h')
       if(idx .eq. 1) then
          write(*,'("============= _qrm_testing usage =============")')
          stop
       end if
    end if

    i = 1
    do
       if(i .gt. nargs) exit

       call get_command_argument(i,value=str,length=len)
       select case(str(1:2))
       case('-t')
          i = i+1
          if(i .gt. nargs) exit
          call get_command_argument(i,value=str,length=len)
          read(str,*)t
       case('-m')
          i = i+1
          if(i .gt. nargs) exit
          call get_command_argument(i,value=str,length=len)
          read(str,*)m
       case('-c')
          i = i+1
          if(i .gt. nargs) exit
          call get_command_argument(i,value=str,length=len)
          read(str,*)c
       case default
          write(*,'("Unrecognized option (try with -h)")')
       end select
       i = i+1
    end do

    return

  end subroutine get_args
  


end program _qrm_testing
