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
!> @file qrm_test_err.F90
!! This file contains tests for the error handling
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################

function _qrm_test_err(c) result(ok)
  use _qrm_testing_mod, protect => _qrm_test_err
  implicit none

  integer :: c
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
  
  if(cases(1 )) ok = case1 ()  .and. ok 
  if(cases(2 )) ok = case2 ()  .and. ok 
  if(cases(3 )) ok = case3 ()  .and. ok 
  if(cases(4 )) ok = case4 ()  .and. ok 
  if(cases(5 )) ok = case5 ()  .and. ok 
  if(cases(6 )) ok = case6 ()  .and. ok 
  if(cases(7 )) ok = case7 ()  .and. ok 
  if(cases(8 )) ok = case8 ()  .and. ok 
  if(cases(9 )) ok = case9 ()  .and. ok 
  if(cases(10)) ok = case10()  .and. ok 
  if(cases(11)) ok = case11()  .and. ok 
  if(cases(12)) ok = case12()  .and. ok 

  return


contains

  !> @brief Check for error 1: unsupported matrix format
  function case1() result(ok)
    use _qrm_spmat_mod
    implicit none

    logical               :: ok
    
    type(_qrm_spmat_type) :: qrm_spmat
    integer               :: info

    call _qrm_spmat_alloc(qrm_spmat, 1, 1, 1, 'xyz', info)
    call _qrm_prnt_testmesg(1, "error", 1, 1, -1, info.eq.1)
    call _qrm_spmat_destroy(qrm_spmat)
    ok = info.eq.1
    return

  end function case1
  

  !> @brief Check for error 4: allocating an already allocated
  function case2() result(ok)
    use qrm_mem_mod
    implicit none
    logical                :: ok

    _qrm_data, allocatable :: a(:)
    _qrm_data, pointer     :: p(:)
    integer :: infoa, infop

    nullify(p)
    
    infoa = 0
    infop = 0
    call qrm_alloc(a, 10)
    call qrm_alloc(a, 10, infoa)
    call qrm_alloc(p, 10)
    call qrm_alloc(p, 10, infop)
    call _qrm_prnt_testmesg(1, "error", 2, 1, -1, (infoa .eq. 4).and.(infop.eq.4))
    call qrm_dealloc(a)
    call qrm_dealloc(p)
    ok = (infoa .eq. 4).and.(infop.eq.4)

    return

  end function case2
  


  !> @brief Check for error 8: invalid or unprovided colperm
  function case3() result(ok)
    use _qrm_mod
    implicit none
    logical                :: ok
    
    type(_qrm_spmat_type), pointer :: qrm_spmat

    integer :: info, i

    qrm_spmat => _qrm_get_test_mat(1)
    call qrm_set(qrm_spmat,'qrm_ordering',qrm_given_)
    ! subcase one, cperm_in not provided
    call qrm_analyse(qrm_spmat, info=info)
    ok = info .eq. 8
    call _qrm_prnt_testmesg(1, "error", 3, 1, 1, ok)

    ! subcase two, cperm_in is not valid
    call qrm_alloc(qrm_spmat%cperm_in, qrm_spmat%n)
    qrm_spmat%cperm_in = (/(i,i=0,qrm_spmat%n-1)/)
    call qrm_analyse(qrm_spmat, info=info)
    call _qrm_prnt_testmesg(1, "error", 3, 2, 1, info .eq. 8)
    call qrm_dealloc(qrm_spmat%cperm_in)
    nullify(qrm_spmat)
    ok = info.eq.8
    return

  end function case3
  

  !> @brief Check for error 9: unknown ordering method
  function case4() result(ok)
    use _qrm_mod
    implicit none
    logical                :: ok
    
    type(_qrm_spmat_type), pointer :: qrm_spmat

    integer :: info

    qrm_spmat => _qrm_get_test_mat(1)
    call qrm_set(qrm_spmat,'qrm_ordering',10)

    call qrm_analyse(qrm_spmat, info=info)
    call _qrm_prnt_testmesg(1, "error", 4, 1, 1, info .eq. 9)
    ok = info.eq.9
    nullify(qrm_spmat)
    return
  end function case4



  !> @brief Check for error 13: facto before analysis
  function case5() result(ok)
    use _qrm_mod
    implicit none
    logical                :: ok
    
    type(_qrm_spmat_type), pointer :: qrm_spmat

    integer :: info

    qrm_spmat => _qrm_get_test_mat(1)
    call qrm_factorize(qrm_spmat, info=info)
    call _qrm_prnt_testmesg(1, "error", 5, 1, 1, info .eq. 13)
    ok = info.eq.13
    nullify(qrm_spmat)
    return
  end function case5


  !> @brief Check for error 14: solve before facto
  function case6() result(ok)
    use _qrm_mod
    implicit none
    logical                :: ok
    
    type(_qrm_spmat_type), pointer :: qrm_spmat
    _qrm_data, allocatable :: b(:), x(:)
    integer   :: info

    qrm_spmat => _qrm_get_test_mat(1)
    call qrm_alloc(b, qrm_spmat%m)
    call qrm_alloc(x, qrm_spmat%n)

    call qrm_apply(qrm_spmat, 'n', b, info)
    call _qrm_prnt_testmesg(1, "error", 6, 1, 1, info .eq. 14)
    ok = info.eq.14
    
    call qrm_solve(qrm_spmat, 'n', b, x, info)
    call _qrm_prnt_testmesg(1, "error", 6, 2, 1, info .eq. 14)
    ok = ok .and. (info.eq.14)
    nullify(qrm_spmat)
    call qrm_dealloc(b)
    call qrm_dealloc(x)
    return
  end function case6


  !> @brief Check for error 15: norm not implemented
  function case7() result(ok)
    use _qrm_mod
    implicit none
    logical                :: ok
    
    _qrm_data, allocatable :: b(:), x(:)
    integer   :: info
    _qrm_real :: nrm

    call qrm_alloc(x, 10)
    call _qrm_vecnrm(x, 10, 'a', nrm, info)
    call _qrm_prnt_testmesg(1, "error", 7, 1, -1, info .eq. 15)
    call qrm_dealloc(x)
    ok = info.eq.15
    return
  end function case7



  !> @brief Check for error 16: unknown ordering method
  function case8() result(ok)
    use _qrm_mod
    implicit none
    logical                :: ok
    
    type(_qrm_spmat_type), pointer :: qrm_spmat

    integer :: info

    qrm_spmat => _qrm_get_test_mat(1)

    call qrm_set(qrm_spmat,'qrm_ordering',qrm_colamd_)
    call qrm_analyse(qrm_spmat, info=info)
#if defined(have_colamd)
    call _qrm_prnt_testmesg(1, "error", 8, 1, 1, info .eq. 0)
    ok = info.eq.0
#else
    call _qrm_prnt_testmesg(1, "error", 8, 1, 1, info .eq. 16)
    ok = info.eq.16
#endif

    call qrm_set(qrm_spmat,'qrm_ordering',qrm_metis_)
    call qrm_analyse(qrm_spmat, info=info)
#if defined(have_metis)
    call _qrm_prnt_testmesg(1, "error", 8, 2, 1, info .eq. 0)
    ok = ok .and. (info.eq.0)
#else
    call _qrm_prnt_testmesg(1, "error", 8, 2, 1, info .eq. 16)
    ok = ok .and. (info.eq.16)
#endif

    call qrm_set(qrm_spmat,'qrm_ordering',qrm_scotch_)
    call qrm_analyse(qrm_spmat, info=info)
#if defined(have_scotch)
    call _qrm_prnt_testmesg(1, "error", 8, 3, 1, info .eq. 0)
    ok = ok .and. (info.eq.0)
#else
    call _qrm_prnt_testmesg(1, "error", 8, 3, 1, info .eq. 16)
    ok = ok .and. (info.eq.16)
#endif

    call qrm_cntl_init(qrm_spmat)
    nullify(qrm_spmat)
    return
  end function case8


  !> @brief Check for error 23: wrong argument to set/get
  function case9() result(ok)
    use _qrm_mod
    implicit none
    logical                :: ok
    
    type(_qrm_spmat_type), pointer :: qrm_spmat

    integer :: info

    qrm_spmat => _qrm_get_test_mat(1)
    call qrm_set(qrm_spmat,'qrm_pippo',qrm_colamd_, info)
    call _qrm_prnt_testmesg(1, "error", 9, 1, 1, info .eq. 23)
    ok = info.eq.23

    nullify(qrm_spmat)
    return
  end function case9


  !> @brief Check for error 28: incorrect blocking sizes
  function case10() result(ok)
    use _qrm_mod
    implicit none
    logical                :: ok
    
    type(_qrm_spmat_type), pointer :: qrm_spmat

    integer                        :: info

    qrm_spmat => _qrm_get_test_mat(1)

    call qrm_set(qrm_spmat,'qrm_nb',-10)
    call qrm_set(qrm_spmat,'qrm_ib',-15)
    call qrm_analyse(qrm_spmat, info=info)
    call _qrm_prnt_testmesg(1, "error", 10, 1, 1, info .eq. 28)
    call qrm_adata_destroy(qrm_spmat%adata)
    ok = info.eq.28

    call qrm_set(qrm_spmat,'qrm_mb',15)
    call qrm_set(qrm_spmat,'qrm_nb',10)
    call qrm_analyse(qrm_spmat, info=info)
    call _qrm_prnt_testmesg(1, "error", 10, 2, 1, info .eq. 28)
    call qrm_adata_destroy(qrm_spmat%adata)
    ok = ok .and. (info.eq.28)

    call qrm_set(qrm_spmat,'qrm_nb',20)
    call qrm_set(qrm_spmat,'qrm_mb',10)
    call qrm_analyse(qrm_spmat, info=info)
    call _qrm_prnt_testmesg(1, "error", 10, 3, 1, info .eq. 28)
    call qrm_adata_destroy(qrm_spmat%adata)
    ok = ok .and. (info.eq.28)

    call qrm_cntl_init(qrm_spmat)
    nullify(qrm_spmat)
    return
  end function case10



  !> @brief Check for error 29: wrong m/n/nz
  function case11() result(ok)
    use _qrm_mod
    implicit none
    logical                :: ok
    
    type(_qrm_spmat_type), pointer :: qrm_spmat

    integer                        :: info, m, n, nz

    qrm_spmat => _qrm_get_test_mat(1)

    m  = qrm_spmat%m
    n  = qrm_spmat%n
    nz = qrm_spmat%nz

    qrm_spmat%m  = -1

    call qrm_analyse(qrm_spmat, info=info)
    call _qrm_prnt_testmesg(1, "error", 11, 1, 1, info .eq. 29)
    call qrm_adata_destroy(qrm_spmat%adata)
    ok = info.eq.29
    
    qrm_spmat%m  = m
    qrm_spmat%n  = n
    qrm_spmat%nz = m*n+1
    call qrm_analyse(qrm_spmat, info=info)
    call _qrm_prnt_testmesg(1, "error", 11, 2, 1, info .eq. 29)
    call qrm_adata_destroy(qrm_spmat%adata)
    ok = ok .and. (info.eq.29)

    qrm_spmat%m  = m
    qrm_spmat%n  = n
    qrm_spmat%nz = nz

    nullify(qrm_spmat)
    return
  end function case11



  !> @brief Check for error 30: apply when H is discarded and solve when R is discarded
  function case12() result(ok)
    use _qrm_mod
    implicit none
    logical                :: ok
    
    type(_qrm_spmat_type), pointer :: qrm_spmat
    _qrm_data, allocatable         :: b(:), x(:)
    integer                        :: info

    qrm_spmat => _qrm_get_test_mat(1)
    call qrm_alloc(b, qrm_spmat%m)
    call qrm_alloc(x, qrm_spmat%n)

    call qrm_set(qrm_spmat, 'qrm_keeph', qrm_no_)
    call qrm_analyse(qrm_spmat)
    call qrm_factorize(qrm_spmat)
    call qrm_apply(qrm_spmat, _qrm_transp, b, info=info)
    call _qrm_prnt_testmesg(1, "error", 12, 1, 1, info .eq. 30)
    ok = info.eq.30

    call _qrm_spmat_cleanup(qrm_spmat)
    call qrm_cntl_init(qrm_spmat)
    call qrm_set(qrm_spmat, 'qrm_keeph', -1)
    call qrm_analyse(qrm_spmat)
    call qrm_factorize(qrm_spmat)
    call qrm_solve(qrm_spmat, 'n', b, x, info=info)
    call _qrm_prnt_testmesg(1, "error", 12, 1, 1, info .eq. 30)
    ok = ok .and. (info.eq.30)

    call _qrm_spmat_cleanup(qrm_spmat)
    call qrm_cntl_init(qrm_spmat)
    nullify(qrm_spmat)
    call qrm_dealloc(b)
    call qrm_dealloc(x)
    return
  end function case12


end function _qrm_test_err
