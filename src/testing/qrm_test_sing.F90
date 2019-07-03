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
!> @file qrm_coverage.F90
!! This file contains coverage tests for the singletons detection
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################

subroutine _qrm_test_sing(c, m)
  use _qrm_coverage_mod, protect => _qrm_test_sing
  implicit none

  integer :: c, m

  integer, parameter :: ncases=15
  logical :: cases(ncases)

  if(c .eq. -1) then
     cases = .true.
  else if (c .le. ncases) then
     cases = .false.
     cases(c) = .true.
  end if

  if(cases(1 )) call case1 (m)

  return


contains

  !> @brief Test the singletons feature
  subroutine case1(m)
    use _qrm_mod
    implicit none
    
    integer :: m

    type(_qrm_spmat_type), pointer :: qrm_spmat
    integer                        :: info
    _qrm_data, pointer             :: x(:), b(:), r(:)
    _qrm_real                      :: err
    integer                        :: i

    if((m.ne.-1) .and. &
         & (m.lt.11) .and. &
         & (m.gt.size(matrices))) then
       write(*,'("Matrix ",i2," is not available for this test")')m
       return
    end if

    ! loop over all the file-matrices
    do i=11, size(matrices)
       if((m.ne.-1) .and. &
            & (m.ne.i)) cycle
       
       ! get the data
       qrm_spmat => _qrm_get_test_mat(i)
       b         => _qrm_get_test_b(i)
       x         => _qrm_get_test_x(i)
       r         => _qrm_get_test_r(i)
       
       ! set the singletons option
       call qrm_set(qrm_spmat,'qrm_sing',qrm_yes_)
       
       ! solve and get back the error
       err = _qrm_test_solve(qrm_spmat, b, x, r)
       
       ! check if anything went wrong
       call qrm_err_get(info)
       call qrm_flush_err_stack(.false.)
       
       ! print message
       call _qrm_prnt_testmesg(4, "sing.", 1, 1, i, (info .eq. 0) .and. (err .lt. eps))
       
       ! put cntl back in its initial state
       call _qrm_cntl_init(qrm_spmat)
    end do
    
    return

  end subroutine case1
  



end subroutine _qrm_test_sing
