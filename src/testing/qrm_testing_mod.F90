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
!> @file qrm_testing_mod.F90
!! Various tools for the testing test
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################

module _qrm_testing_mod
  use _qrm_spmat_mod

  interface 
     function _qrm_test_err(c) result(ok)
       integer :: c
       logical :: ok
     end function _qrm_test_err
  end interface

  interface 
     function _qrm_test_ord(c, m) result(ok)
       integer :: c, m
       logical :: ok
     end function _qrm_test_ord
  end interface

  interface 
     function _qrm_test_c(c, m) result(ok)
       integer :: c, m
       logical :: ok
     end function _qrm_test_c
  end interface

  interface 
     function _qrm_test_multiples(c, m) result(ok)
       integer :: c, m
       logical :: ok
     end function _qrm_test_multiples
  end interface

  interface 
     function _qrm_test_facto(c, m) result(ok)
       integer :: c, m
       logical :: ok
     end function _qrm_test_facto
  end interface

  interface 
     function _qrm_test_methods(c, m) result(ok)
       integer :: c, m
       logical :: ok
     end function _qrm_test_methods
  end interface

  interface 
     function _qrm_test_pipe(c, m) result(ok)
       integer :: c, m
       logical :: ok
     end function _qrm_test_pipe
  end interface

  type tmat
     character(len=50)                  :: mfile
     type(_qrm_spmat_type), allocatable :: qrm_spmat
  end type tmat

  type(tmat), target, allocatable :: matrices(:)
  integer :: iseed(4)

#if defined(dprec) || defined(zprec)
  _qrm_real, parameter :: eps=10d-9
#elif defined(sprec) || defined(cprec)
  _qrm_real, parameter :: eps=10e-3
#endif


contains

  subroutine _qrm_prnt_testmesg(t, name, c, s, m, res)
    
    integer          :: t, m, c, s
    character(len=5) :: name
    logical          :: res

    if(res) then
       write(*,'("TEST ",i2," (",a5,")     CASE ",i2,"     SUBCASE ",i2,"     MATRIX ",i2," : OK")')t,name,c,s,m
    else
       write(*,'("TEST ",i2," (",a5,")     CASE ",i2,"     SUBCASE ",i2,"     MATRIX ",i2," : FAILED")')t,name,c,s,m
    end if

    return
  end subroutine _qrm_prnt_testmesg
  


  subroutine read_matfile()
    use _qrm_spmat_mod
    implicit none

    integer :: info, nmats, i

    open(4,file='matfile.txt', status='OLD', action='READ', iostat=info)
    read(4,*)nmats
    
    allocate(matrices(nmats+10))

    do i=1, nmats
       read(4,'(a50)')matrices(10+i)%mfile
    end do

    iseed = (/1,1,1,1/)
    return

  end subroutine read_matfile
  

  function _qrm_get_test_mat(m)
    use _qrm_utils_mod
    implicit none
    integer :: m
    type(_qrm_spmat_type), pointer :: _qrm_get_test_mat
    integer :: nm

    select case (m)
    case(1)
       if(.not. allocated(matrices(1)%qrm_spmat)) then
          allocate(matrices(1)%qrm_spmat)
          call _qrm_spmat_init(matrices(1)%qrm_spmat)
          call _qrm_spmat_alloc(matrices(1)%qrm_spmat, 5, 4, 4, 'coo')
          matrices(1)%qrm_spmat%irn=(/1, 2, 3, 4, 4/)
          matrices(1)%qrm_spmat%jcn=(/1, 2, 3, 4, 3/)
          call _xlarnv(2,iseed,matrices(1)%qrm_spmat%nz, matrices(1)%qrm_spmat%val)
       end if
       _qrm_get_test_mat => matrices(1)%qrm_spmat
    case(11:)
       if(m.gt.size(matrices)) then
          write(*,'("Matrix ",i2," does not exist")')m
          return
       end if
       if(.not. allocated(matrices(m)%qrm_spmat)) then
          allocate(matrices(m)%qrm_spmat)
          call _qrm_spmat_init(matrices(m)%qrm_spmat)
          call qrm_readmat(matrices(m)%mfile, matrices(m)%qrm_spmat, .true.)
          ! else
          ! call _qrm_spmat_destroy(matrices(m)%qrm_spmat)
       end if
       ! call _qrm_spmat_cleanup(matrices(m)%qrm_spmat)
       ! call _qrm_cntl_init(matrices(m)%qrm_spmat)
       _qrm_get_test_mat => matrices(m)%qrm_spmat
    case default
       write(*,'("Matrix ",i2," does not exist")')m
       return
    end select

    return

  end function _qrm_get_test_mat

end module _qrm_testing_mod
