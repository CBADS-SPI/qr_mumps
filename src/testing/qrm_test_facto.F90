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
!> @file qrm_test_facto.F90
!! This file contains coverage tests for the factorization
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################

function _qrm_test_facto(c, m) result(ok)
  use _qrm_testing_mod, protect => _qrm_test_facto
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
  if(cases(3 )) ok = case3 (m) .and. ok 

  return


contains

  !> @brief Test the factorization 
  function case1(m) result(ok)
    use _qrm_mod
    implicit none
    
    integer :: m
    logical :: ok
    
    type(_qrm_spmat_type), pointer :: qrm_spmat
    integer                        :: info
    _qrm_data, allocatable         :: x(:), b(:), r(:)
    _qrm_real                      :: err, anrm, rnrm, bnrm, xnrm, onrm
    integer                        :: i, j
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
       if(info.eq.0) call qrm_alloc(b, qrm_spmat%m, info)
       if(info.eq.0) call qrm_alloc(r, qrm_spmat%m, info)
       if(info.eq.0) call qrm_alloc(x, qrm_spmat%n, info)
       if(info.ne.0) goto 9999
       call _xlarnv(2, iseed, size(b), b(1))
       r = b
       
       if(qrm_spmat%m .ge. qrm_spmat%n) then
          transp='n'
       else
          transp=_qrm_transp
       end if
       
       call _xlarnv(2, iseed, size(b), b(1))
       r = b
       
       call qrm_analyse(qrm_spmat, transp, info)
       if(info.ne.0) goto 9999
       call qrm_factorize(qrm_spmat, transp, info)
       if(info.ne.0) goto 9999
       
       if(transp .eq. 'n') then
          call qrm_apply(qrm_spmat, _qrm_transp, b, info)
          if(info.ne.0) goto 9999
          call qrm_solve(qrm_spmat, 'n', b, x, info)
          if(info.ne.0) goto 9999
       else if(transp .eq. _qrm_transp) then
          call qrm_solve(qrm_spmat, _qrm_transp, b, x, info)
          if(info.ne.0) goto 9999
          call qrm_apply(qrm_spmat, 'n', x, info)
          if(info.ne.0) goto 9999
       end if
       
       call qrm_residual_norm(qrm_spmat, r, x, rnrm)
       call qrm_vecnrm(x, size(x,1), '2', xnrm)
       call qrm_vecnrm(b, size(b,1), '2', bnrm)
       call qrm_matnrm(qrm_spmat, 'f', anrm)
       call qrm_residual_orth(qrm_spmat, r, onrm)   
       
       if(transp .eq. 'n') then
          if(rnrm.lt.eps) then
             err = rnrm/anrm
          else
             err = onrm/rnrm
          end if
       else if(transp .eq. _qrm_transp) then
          err = rnrm/anrm
       end if
       
       ! print message
9999   continue
       call _qrm_prnt_testmesg(5, "facto", 1, 1, i, (info .eq. 0) .and. (err .lt. eps))
       ok = ok .and. (info .eq. 0) .and. (err .lt. eps)
       
       ! put cntl back in its initial state
       call _qrm_spmat_cleanup(qrm_spmat)
       call _qrm_cntl_init(qrm_spmat)

       call qrm_dealloc(b)
       call qrm_dealloc(r)
       call qrm_dealloc(x)

    end do
    

    return

  end function case1
  

  !> @brief Test the factorization with different
  !! mb, nb, ib and bh parameters
  function case2(m) result(ok)
    use _qrm_mod
    implicit none
    
    integer :: m
    logical :: ok

    type(_qrm_spmat_type), pointer :: qrm_spmat
    integer                        :: info
    _qrm_data, allocatable         :: x(:), b(:), r(:)
    _qrm_real                      :: err, anrm, rnrm, bnrm, xnrm, onrm
    integer                        :: i, j, p
    character                      :: transp
    integer, parameter             :: np=6
    integer                        :: params(4,np)

    params = reshape((/  64, 64, 32, 0, &
                      & 128, 64, 32, 0, &
                      &  64, 64, 32, 2, &
                      & -64, 64, 64, 0, &
                      & 160, 80, 40, 1, &
                      & 512, 64, 32, 0 /), shape(params))
    
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
       
       ! get the data
       qrm_spmat => _qrm_get_test_mat(i)
       if(qrm_spmat%m .ge. qrm_spmat%n) then
          transp='n'
       else
          transp=_qrm_transp
       end if

       do p=1, 6
          info = 0

          if(info.eq.0) call qrm_alloc(b, qrm_spmat%m, info)
          if(info.eq.0) call qrm_alloc(r, qrm_spmat%m, info)
          if(info.eq.0) call qrm_alloc(x, qrm_spmat%n, info)
          if(info.ne.0) goto 9999
          call _xlarnv(2, iseed, size(b), b(1))
          r = b

          call qrm_set(qrm_spmat, 'qrm_mb_', params(1,p))
          call qrm_set(qrm_spmat, 'qrm_nb_', params(2,p))
          call qrm_set(qrm_spmat, 'qrm_ib_', params(3,p))
          call qrm_set(qrm_spmat, 'qrm_bh_', params(4,p))
          
          call _xlarnv(2, iseed, size(b), b(1))
          r = b
          
          call qrm_analyse(qrm_spmat, transp, info)
          if(info.ne.0) goto 9999
          call qrm_factorize(qrm_spmat, transp, info)
          if(info.ne.0) goto 9999

          if(transp .eq. 'n') then
             call qrm_apply(qrm_spmat, _qrm_transp, b, info)
             if(info.ne.0) goto 9999
             call qrm_solve(qrm_spmat, 'n', b, x, info)
          if(info.ne.0) goto 9999
          else if(transp .eq. _qrm_transp) then
             call qrm_solve(qrm_spmat, _qrm_transp, b, x, info)
             if(info.ne.0) goto 9999
             call qrm_apply(qrm_spmat, 'n', x, info)
             if(info.ne.0) goto 9999
          end if
          
          call qrm_residual_norm(qrm_spmat, r, x, rnrm)
          call qrm_vecnrm(x, size(x,1), '2', xnrm)
          call qrm_vecnrm(b, size(b,1), '2', bnrm)
          call qrm_matnrm(qrm_spmat, 'f', anrm)
          call qrm_residual_orth(qrm_spmat, r, onrm)   
          
          if(transp .eq. 'n') then
             if(rnrm.lt.eps) then
                err = rnrm/anrm
             else
                err = onrm/rnrm
             end if
          else if(transp .eq. _qrm_transp) then
             err = rnrm/anrm
          end if
          
9999      continue
          ! print message
          call _qrm_prnt_testmesg(5, "facto", 2, p, i, (info .eq. 0) .and. (err .lt. eps))
          ok = ok .and. (info .eq. 0) .and. (err .lt. eps)
          
          ! cleanup the matrix and put cntl back in its initial state
          call _qrm_spmat_cleanup(qrm_spmat)
          call _qrm_cntl_init(qrm_spmat)
          call qrm_dealloc(b)
          call qrm_dealloc(r)
          call qrm_dealloc(x)

       end do
    end do
    return

  end function case2


  !> @brief Test the factorization with different
  !! memory bounds
  function case3(m) result(ok)
    use _qrm_mod
    implicit none
    
    integer :: m
    logical :: ok

    type(_qrm_spmat_type), pointer :: qrm_spmat
    integer                        :: info
    _qrm_data, allocatable         :: x(:), b(:), r(:)
    _qrm_real                      :: err, anrm, rnrm, bnrm, xnrm, onrm
    integer                        :: i, j, p
    character                      :: transp
    real(kind(1.d0))               :: memrelax


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
       
       ! get the data
       qrm_spmat => _qrm_get_test_mat(i)
       if(qrm_spmat%m .ge. qrm_spmat%n) then
          transp='n'
       else
          transp=_qrm_transp
       end if

       do p=1, 3
          memrelax = 1.d0+(p-1)*0.25d0
          info = 0

          if(info.eq.0) call qrm_alloc(b, qrm_spmat%m, info)
          if(info.eq.0) call qrm_alloc(r, qrm_spmat%m, info)
          if(info.eq.0) call qrm_alloc(x, qrm_spmat%n, info)
          if(info.ne.0) goto 9999
          call _xlarnv(2, iseed, size(b), b(1))
          r = b

          call qrm_set(qrm_spmat, 'qrm_mem_relax_', memrelax)
          
          call _xlarnv(2, iseed, size(b), b(1))
          r = b
          
          call qrm_analyse(qrm_spmat, transp, info)
          if(info.ne.0) goto 9999
          call qrm_factorize(qrm_spmat, transp, info)
          if(info.ne.0) goto 9999

          if(transp .eq. 'n') then
             call qrm_apply(qrm_spmat, _qrm_transp, b, info)
             if(info.ne.0) goto 9999
             call qrm_solve(qrm_spmat, 'n', b, x, info)
          if(info.ne.0) goto 9999
          else if(transp .eq. _qrm_transp) then
             call qrm_solve(qrm_spmat, _qrm_transp, b, x, info)
             if(info.ne.0) goto 9999
             call qrm_apply(qrm_spmat, 'n', x, info)
             if(info.ne.0) goto 9999
          end if
          
          call qrm_residual_norm(qrm_spmat, r, x, rnrm)
          call qrm_vecnrm(x, size(x,1), '2', xnrm)
          call qrm_vecnrm(b, size(b,1), '2', bnrm)
          call qrm_matnrm(qrm_spmat, 'f', anrm)
          call qrm_residual_orth(qrm_spmat, r, onrm)   
          
          if(transp .eq. 'n') then
             if(rnrm.lt.eps) then
                err = rnrm/anrm
             else
                err = onrm/rnrm
             end if
          else if(transp .eq. _qrm_transp) then
             err = rnrm/anrm
          end if
          
9999      continue
          ! print message
          call _qrm_prnt_testmesg(5, "facto", 3, p, i, (info .eq. 0) .and. (err .lt. eps))
          ok = ok .and. (info .eq. 0) .and. (err .lt. eps)
          
          ! cleanup the matrix and put cntl back in its initial state
          call _qrm_spmat_cleanup(qrm_spmat)
          call _qrm_cntl_init(qrm_spmat)
          call qrm_dealloc(b)
          call qrm_dealloc(r)
          call qrm_dealloc(x)

       end do
    end do
    return

  end function case3

  

end function _qrm_test_facto
