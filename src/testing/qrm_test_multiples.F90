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
!> @file qrm_test_multiples.F90
!! This file contains coverage tests for multiple consecutive call of the same routine on the same data
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################

function _qrm_test_multiples(c, m) result(ok)
  use _qrm_testing_mod, protect => _qrm_test_multiples
  implicit none

  integer :: c, m
  logical :: ok
  
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

  !> @brief Test multiple call to the analysis routine
  function case1(m) result(ok)
    use _qrm_mod
    implicit none
    
    integer :: m
    logical :: ok
    
    type(_qrm_spmat_type), pointer :: qrm_spmat
    integer                        :: info
    integer                        :: i
    character                      :: transp 

    
    if((m.ne.-1) .and. &
         & (m.lt.11) .and. &
         & (m.gt.size(matrices))) then
       write(*,'("Matrix ",i2," is not available for this test")')m
       return
    end if

    ok = .true.
    
    ! loop over all the file-matrices
    do i=11, size(matrices)
       if((m.ne.-1) .and. &
            & (m.ne.i)) cycle
       
       info = 0

       ! get the data
       qrm_spmat => _qrm_get_test_mat(i)
       
       if(qrm_spmat%m .ge. qrm_spmat%n) then
          transp='n'
       else
          transp=_qrm_transp
       end if
       
       call qrm_analyse(qrm_spmat, transp, info)
       if(info.ne.0) goto 9999
       
       call qrm_analyse(qrm_spmat, transp, info)
       if(info.ne.0) goto 9999
       
       ! print message
9999   continue
       call _qrm_prnt_testmesg(4, "mult.", 1, 1, i, info .eq. 0)
       ok = ok .and. (info.eq.0)
       
       ! put cntl back in its initial state
       call _qrm_spmat_cleanup(qrm_spmat)
       call _qrm_cntl_init(qrm_spmat)

    end do
    

    return

  end function case1



  !> @brief Test multiple call to the factorizatopn routine
  function case2(m) result(ok)
    use _qrm_mod
    implicit none
    
    integer :: m
    logical :: ok

    type(_qrm_spmat_type), pointer :: qrm_spmat
    integer                        :: info
    integer                        :: i
    character                      :: transp 

    
    if((m.ne.-1) .and. &
         & (m.lt.11) .and. &
         & (m.gt.size(matrices))) then
       write(*,'("Matrix ",i2," is not available for this test")')m
       return
    end if

    ok = .true.
    
    ! loop over all the file-matrices
    do i=11, size(matrices)
       if((m.ne.-1) .and. &
            & (m.ne.i)) cycle
       
       info = 0

       ! get the data
       qrm_spmat => _qrm_get_test_mat(i)
       
       if(qrm_spmat%m .ge. qrm_spmat%n) then
          transp='n'
       else
          transp=_qrm_transp
       end if
       
       call qrm_analyse(qrm_spmat, transp, info)
       if(info.ne.0) goto 9999
       
       call qrm_factorize(qrm_spmat, transp, info)
       if(info.ne.0) goto 9999
       
       call qrm_factorize(qrm_spmat, transp, info)
       if(info.ne.0) goto 9999
       
       ! print message
9999   continue
       call _qrm_prnt_testmesg(4, "mult.", 2, 1, i, info .eq. 0)
       ok = ok .and. (info.eq.0)
       
       ! put cntl back in its initial state
       call _qrm_spmat_cleanup(qrm_spmat)
       call _qrm_cntl_init(qrm_spmat)

    end do
    

    return

  end function case2
  

  

end function _qrm_test_multiples
