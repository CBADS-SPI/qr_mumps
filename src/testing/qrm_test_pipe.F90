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
!> @file qrm_test_pipe.F90
!! This file contains coverage tests for the factorization
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################

function _qrm_test_pipe(c, m) result(ok)
  use _qrm_testing_mod, protect => _qrm_test_pipe
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

  !> @brief Test the factorization with different
  !! numbers of threads
  function case1(m) result(ok)
    use _qrm_mod
    implicit none

    integer :: m
    logical :: ok

    type(_qrm_spmat_type), pointer :: qrm_spmat
    type(qrm_dscr_type)            :: qrm_dscr
    integer                        :: info
    _qrm_data, allocatable         :: x(:), b(:), r(:)
    type(_qrm_rhs_type)            :: x_rhs, b_rhs
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

       call _qrm_rhs_init(b_rhs, b)
       call _qrm_rhs_init(x_rhs, x)

       if(qrm_spmat%m .ge. qrm_spmat%n) then
          transp='n'
       else
          transp=_qrm_transp
       end if

       call qrm_dscr_init(qrm_dscr)

       call qrm_analyse_async(qrm_dscr, qrm_spmat, transp, info)
       if(info.ne.0) goto 9998
       call qrm_factorize_async(qrm_dscr, qrm_spmat, transp, info)
       if(info.ne.0) goto 9998

       if(transp .eq. 'n') then
          call qrm_apply_async(qrm_dscr, qrm_spmat, _qrm_transp, b_rhs, info)
          if(info.ne.0) goto 9998
          call qrm_solve_async(qrm_dscr, qrm_spmat, 'n', b_rhs, x_rhs, info)
          if(info.ne.0) goto 9998
       else if(transp .eq. _qrm_transp) then
          call qrm_solve_async(qrm_dscr, qrm_spmat, _qrm_transp, b_rhs, x_rhs, info)
          if(info.ne.0) goto 9998
          call qrm_apply_async(qrm_dscr, qrm_spmat, 'n', x_rhs, info)
          if(info.ne.0) goto 9998
       end if

9998   continue
       call qrm_barrier(qrm_dscr)
       info = merge(info, qrm_dscr%err_status, info.ne.0)
       if(info.ne.0) goto 9999


       call qrm_dscr_destroy(qrm_dscr)
       call _qrm_rhs_destroy(b_rhs)
       call _qrm_rhs_destroy(x_rhs)

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

9999   continue
       ! print message
       call _qrm_prnt_testmesg(7, "pipe.", 1, 1, i, (info .eq. 0) .and. (err .lt. eps))
       ok = ok .and. (info .eq. 0) .and. (err .lt. eps)

       ! cleanup the matrix and put cntl back in its initial state
       call _qrm_spmat_cleanup(qrm_spmat)
       call _qrm_cntl_init(qrm_spmat)

       call qrm_dealloc(b)
       call qrm_dealloc(r)
       call qrm_dealloc(x)

    end do

    return

  end function case1

  !> @brief Test the factorization with different
  !! numbers of threads
  function case2(m) result(ok)
    use _qrm_mod
    implicit none

    integer :: m
    logical :: ok

    type(qrm_dscr_type)            :: qrm_dscr1, qrm_dscr2
    type(_qrm_spmat_type), pointer :: qrm_spmat1
    type(_qrm_spmat_type)          :: qrm_spmat2
    integer                        :: info
    _qrm_data, allocatable         :: x1(:), b1(:), r1(:), x2(:), b2(:), r2(:)
    type(_qrm_rhs_type)            :: x1_rhs, b1_rhs, x2_rhs, b2_rhs
    _qrm_real                      :: err, anrm, rnrm, bnrm, xnrm, onrm
    integer                        :: i, j
    character                      :: transp

    if((m.ne.-1) .and. &
         & (m.lt.11) .and. &
         & (m.gt.size(matrices))) then
       write(*,'("Matrix ",i2," is not available for this test")')m
       return
    end if
    err = huge(_qrm_rone)

    ok = .true.
    
    ! loop over all the file-matrices
    do i=11, size(matrices)
       if((m.ne.-1) .and. &
            & (m.ne.i)) cycle

       info = 0

       ! get the data
       qrm_spmat1 => _qrm_get_test_mat(i)
       if(info.eq.0) call qrm_spmat_copy(qrm_spmat1, qrm_spmat2, info=info)
       if(info.eq.0) call qrm_alloc(b1, qrm_spmat1%m, info)
       if(info.eq.0) call qrm_alloc(r1, qrm_spmat1%m, info)
       if(info.eq.0) call qrm_alloc(x1, qrm_spmat1%n, info)
       if(info.ne.0) goto 9999
       call _xlarnv(2, iseed, size(b1), b1(1))
       r1 = b1

       call _qrm_rhs_init(b1_rhs, b1)
       call _qrm_rhs_init(x1_rhs, x1)

       if(info.eq.0) call qrm_alloc(b2, qrm_spmat2%m, info)
       if(info.eq.0) call qrm_alloc(r2, qrm_spmat2%m, info)
       if(info.eq.0) call qrm_alloc(x2, qrm_spmat2%n, info)
       if(info.ne.0) goto 9999
       call _xlarnv(2, iseed, size(b2), b2(1))
       r2 = b2

       call _qrm_rhs_init(b2_rhs, b2)
       call _qrm_rhs_init(x2_rhs, x2)
       call qrm_dscr_init(qrm_dscr1)
       call qrm_dscr_init(qrm_dscr2)

       if(qrm_spmat1%m .ge. qrm_spmat1%n) then
          transp='n'
       else
          transp=_qrm_transp
       end if

       call qrm_analyse_async(qrm_dscr1, qrm_spmat1, transp, info)
       if(info.ne.0) goto 9998
       call qrm_analyse_async(qrm_dscr2, qrm_spmat2, transp, info)
       if(info.ne.0) goto 9998

       call qrm_factorize_async(qrm_dscr1, qrm_spmat1, transp, info)
       if(info.ne.0) goto 9998
       call qrm_factorize_async(qrm_dscr2, qrm_spmat2, transp, info)
       if(info.ne.0) goto 9998

       if(transp .eq. 'n') then
          call qrm_apply_async(qrm_dscr1, qrm_spmat1, _qrm_transp, b1_rhs, info)
          if(info.ne.0) goto 9998
          call qrm_solve_async(qrm_dscr1, qrm_spmat1, 'n', b1_rhs, x1_rhs, info)
          if(info.ne.0) goto 9998
          call qrm_apply_async(qrm_dscr2, qrm_spmat2, _qrm_transp, b2_rhs, info)
          if(info.ne.0) goto 9998
          call qrm_solve_async(qrm_dscr2, qrm_spmat2, 'n', b2_rhs, x2_rhs, info)
          if(info.ne.0) goto 9998
       else if(transp .eq. _qrm_transp) then
          call qrm_solve_async(qrm_dscr1, qrm_spmat1, _qrm_transp, b1_rhs, x1_rhs, info)
          if(info.ne.0) goto 9998
          call qrm_apply_async(qrm_dscr1, qrm_spmat1, 'n', x1_rhs, info)
          if(info.ne.0) goto 9998
          call qrm_solve_async(qrm_dscr2, qrm_spmat2, _qrm_transp, b2_rhs, x2_rhs, info)
          if(info.ne.0) goto 9998
          call qrm_apply_async(qrm_dscr2, qrm_spmat2, 'n', x2_rhs, info)
          if(info.ne.0) goto 9998
       end if

9998   continue
       call qrm_barrier(qrm_dscr1)
       call qrm_barrier(qrm_dscr2)
       info = merge(info, qrm_dscr1%err_status, info.ne.0)
       info = merge(info, qrm_dscr2%err_status, info.ne.0)
       if(info.ne.0) goto 9999

       call qrm_dscr_destroy(qrm_dscr1)
       call qrm_dscr_destroy(qrm_dscr2)

       call _qrm_rhs_destroy(b1_rhs)
       call _qrm_rhs_destroy(x1_rhs)

       call _qrm_rhs_destroy(b2_rhs)
       call _qrm_rhs_destroy(x2_rhs)

       call qrm_residual_norm(qrm_spmat1, r1, x1, rnrm)
       call qrm_vecnrm(x1, size(x1,1), '2', xnrm)
       call qrm_vecnrm(b1, size(b1,1), '2', bnrm)
       call qrm_matnrm(qrm_spmat1, 'f', anrm)
       call qrm_residual_orth(qrm_spmat1, r1, onrm)

       if(transp .eq. 'n') then
          if(rnrm.lt.eps) then
             err = rnrm/anrm
          else
             err = onrm/rnrm
          end if
       else if(transp .eq. _qrm_transp) then
             err = rnrm/anrm
       end if

       call qrm_residual_norm(qrm_spmat2, r2, x2, rnrm)
       call qrm_vecnrm(x2, size(x2,1), '2', xnrm)
       call qrm_vecnrm(b2, size(b2,1), '2', bnrm)
       call qrm_residual_orth(qrm_spmat2, r2, onrm)

       if(transp .eq. 'n') then
          if(rnrm.lt.eps) then
             err = max(rnrm/anrm,err)
          else
             err = max(onrm/rnrm,err)
          end if
       else if(transp .eq. _qrm_transp) then
             err = max(rnrm/anrm,err)
       end if

9999   continue
       ! print message
       call _qrm_prnt_testmesg(7, "pipe.", 2, 1, i, (info .eq. 0) .and. (err .lt. eps))
       ok = ok .and. (info .eq. 0) .and. (err .lt. eps)

       ! cleanup the matrix and put cntl back in its initial state
       call _qrm_spmat_cleanup(qrm_spmat1)
       call _qrm_cntl_init(qrm_spmat1)

       call qrm_spmat_destroy(qrm_spmat2, all=.true.)
       call qrm_dealloc(b1)
       call qrm_dealloc(r1)
       call qrm_dealloc(x1)
       call qrm_dealloc(b2)
       call qrm_dealloc(r2)
       call qrm_dealloc(x2)

    end do

    return

  end function case2


end function _qrm_test_pipe
