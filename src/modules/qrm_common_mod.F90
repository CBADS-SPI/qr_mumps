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
!> @file qrm_common_mod.F90
!! this module contains generic interfaces for all the untyped routines.
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################

#include "qrm_common.h"

!> @brief This module contains the interfaces of all non-typed routines
module qrm_common_mod
  use qrm_const_mod


  !> @brief Generic interface for the @link ::qrm_print_nsteps_tree @endlink,
  !! @link                                  ::qrm_print_elim_tree @endlink and @link
  !!                                        ::qrm_print_asm_tree @endlink routines
  interface qrm_print_tree
     subroutine qrm_print_nsteps_tree(file, adata, small)
       use qrm_adata_mod
       type(qrm_adata_type)                 :: adata
       character                            :: file*(*)
       logical, optional                    :: small
     end subroutine qrm_print_nsteps_tree
  end interface qrm_print_tree

  !> @brief Generic interface for the @link ::qrm_postorder @endlink routine
  interface qrm_postorder
     subroutine qrm_postorder(parent, n, porder, weight, info)
       integer                              :: n
       integer                              :: parent(:), porder(:)
       integer, optional                    :: weight(:)
       integer, optional                    :: info
     end subroutine qrm_postorder
  end interface qrm_postorder

  !> @brief Generic interface for the @link ::qrm_amalg_tree @endlink routine
  interface
     subroutine qrm_amalg_tree(n, parent, rc, porder, nvar, min_var, fill_thresh, info)
       integer                              :: n, min_var
       integer                              :: parent(:), rc(:), porder(:), nvar(:)
       real(kind(1.d0))                     :: fill_thresh
       integer, optional                    :: info
     end subroutine qrm_amalg_tree
  end interface

  !> @brief Generic interface for the @link ::qrm_compress_data @endlink routine
  interface qrm_compress_data
     subroutine qrm_compress_data(adata, porder, parent, rc, stair, n, info)
       use qrm_mem_mod
       use qrm_adata_mod
       type(qrm_adata_type)                 :: adata
       integer                              :: porder(:), parent(:), rc(:), stair(:)
       integer                              :: n
       integer, optional                    :: info
     end subroutine qrm_compress_data
  end interface qrm_compress_data

  !> @brief Generic interface for the @link ::qrm_reorder_tree @endlink routine
  interface qrm_reorder_tree
     subroutine qrm_reorder_tree(adata, info)
       use qrm_adata_mod
       type(qrm_adata_type)                 :: adata
       integer, optional                    :: info
     end subroutine qrm_reorder_tree
  end interface qrm_reorder_tree

  !> @brief Generic interface for the @link ::qrm_prune_tree @endlink routine
  interface qrm_prune_tree
     subroutine qrm_prune_tree(adata, nth, info)
       use qrm_adata_mod
       type(qrm_adata_type)                 :: adata
       integer                              :: nth
       integer, optional                    :: info
     end subroutine qrm_prune_tree
  end interface qrm_prune_tree

  !> @brief Generic interface for the @link ::qrm_prnt_iarray
  !> @endlink, @link                        ::qrm_prnt_sarray @endlink, @link
  !>                                        ::qrm_prnt_darray @endlink, @link ::qrm_prnt_carray @endlink and
  !> @link                                  ::qrm_prnt_zarray @endlink routines
  interface qrm_prnt_array
     subroutine qrm_prnt_iarray(a, lab, unit)
       integer                              :: a(:)
       character                            :: lab*(*)
       integer, optional                    :: unit
     end subroutine qrm_prnt_iarray
     subroutine qrm_prnt_sarray(a, lab, unit)
       real(kind(1.e0))                     :: a(:)
       character                            :: lab*(*)
       integer, optional                    :: unit
     end subroutine qrm_prnt_sarray
     subroutine qrm_prnt_darray(a, lab, unit)
       real(kind(1.d0))                     :: a(:)
       character                            :: lab*(*)
       integer, optional                    :: unit
     end subroutine qrm_prnt_darray
     subroutine qrm_prnt_d2array(a, lab, unit)
       real(kind(1.d0))                     :: a(:,:)
       character                            :: lab*(*)
       integer, optional                    :: unit
     end subroutine qrm_prnt_d2array
     subroutine qrm_prnt_carray(a, lab, unit)
       complex(kind(1.e0))                  :: a(:)
       character                            :: lab*(*)
       integer, optional                    :: unit
     end subroutine qrm_prnt_carray
     subroutine qrm_prnt_zarray(a, lab, unit)
       complex(kind(1.d0))                  :: a(:)
       character                            :: lab*(*)
       integer, optional                    :: unit
     end subroutine qrm_prnt_zarray
  end interface qrm_prnt_array

  !> @brief Generic interface for the @link ::qrm_swtime @endlink routine
  interface qrm_swtime
     function qrm_swtime() bind(c, name='qrm_swtime')
       use iso_c_binding
       real(c_double)                     :: qrm_swtime
     end function qrm_swtime
  end interface qrm_swtime

  !> @brief Generic interface for the @link ::qrm_uwtime @endlink routine
  interface qrm_uwtime
     function qrm_uwtime() bind(c, name='qrm_uwtime')
       use iso_c_binding
       real(c_double)                     :: qrm_uwtime
     end function qrm_uwtime
  end interface qrm_uwtime

  !> @brief Generic interface for the @link ::qrm_msleep @endlink routine
  interface qrm_msleep
     subroutine qrm_msleep(n) bind(c, name='qrm_msleep')
       use iso_c_binding
       integer(c_int), value                :: n
     end subroutine qrm_msleep
  end interface qrm_msleep


  !> @brief Generic interface for the @link ::qrm_check_cperm @endlink routine
  interface qrm_check_cperm
     subroutine qrm_check_cperm(cperm, n, info)
       integer                              :: cperm(:)
       integer                              :: n, info
     end subroutine qrm_check_cperm
  end interface qrm_check_cperm

  !> @brief Generic interface for the
  !! @link                                  ::qrm_count_realflops @endlink
  !! @link                                  ::qrm_count_pureflops @endlink
  interface qrm_count_flops
     module procedure qrm_count_realflops, qrm_count_pureflops
     ! function qrm_count_realflops(m, n, k, op)
     ! real(kind(1.d0))                     :: qrm_count_realflops
     ! integer                              :: m, k, n
     ! character                            :: op*(*)
     ! end function qrm_count_realflops
     ! function qrm_count_pureflops(stair, n, j, nb)
     ! implicit none
     ! real(kind(1.d0))                     :: qrm_count_pureflops
     ! integer                              :: stair(:)
     ! integer                              :: n, nb, j
     ! end function qrm_count_pureflops
  end interface qrm_count_flops

  interface qrm_set
     module procedure  qrm_gseti
  end interface qrm_set

  interface qrm_get
     module procedure  qrm_ggeti, qrm_ggetii
  end interface qrm_get

contains

  subroutine qrm_gseti(string, ival, info)
    use qrm_string_mod
    use qrm_error_mod
    use qrm_mem_mod
    implicit none

    character(len=*)                        :: string
    integer                                 :: ival
    integer, optional                       :: info
    
    character(len=len(string))              :: istring

    ! error management
    integer                                 :: err
    character(len=*), parameter             :: name='qrm_gseti'

    err = 0

    istring = qrm_str_tolower(string)

    if(index(istring,'qrm_eunit') .eq. 1) then
       qrm_eunit = ival
    else if(index(istring,'qrm_ounit') .eq. 1) then
       qrm_ounit = ival
    else if(index(istring,'qrm_dunit') .eq. 1) then
       qrm_dunit = ival
    else if(index(istring,'qrm_max_mem') .eq. 1) then
       qrm_max_mem = ival
    else if(index(istring,'qrm_tot_mem') .eq. 1) then
       qrm_tot_mem = ival
    else
       err = 23
       call qrm_error_print(err, name, aed=string)
       goto 9999
    end if

9999 continue

    if(present(info)) info = err
    return
  end subroutine qrm_gseti

  subroutine qrm_ggeti(string, ival, info)
    use qrm_string_mod
    use qrm_error_mod
    implicit none

    character(len=*)                        :: string
    integer                                 :: ival
    integer, optional                       :: info
    
    integer(kind=8)                         :: iival

    ! error management
    integer                                 :: err
    character(len=*), parameter             :: name='qrm_ggeti'

    err = 0

    call qrm_ggetii(string, iival, err)
    ival = iival

    if(present(info)) info = err
    return

  end subroutine qrm_ggeti

  subroutine qrm_ggetii(string, iival, info)
    use qrm_string_mod
    use qrm_error_mod
    use qrm_mem_mod
    implicit none

    character(len=*)                        :: string
    integer(kind=8)                         :: iival
    integer, optional                       :: info
    
    character(len=len(string))              :: istring

    ! error management
    integer                                 :: err
    character(len=*), parameter             :: name='qrm_ggeti'

    err = 0

    istring = qrm_str_tolower(string)

    if(index(istring,'qrm_max_mem') .eq. 1) then
       iival = qrm_max_mem
    else if(index(istring,'qrm_tot_mem') .eq. 1) then
       iival = qrm_tot_mem
    else if(index(istring,'qrm_ounit') .eq. 1) then
       iival = qrm_ounit
    else if(index(istring,'qrm_eunit') .eq. 1) then
       iival = qrm_eunit
    else if(index(istring,'qrm_dunit') .eq. 1) then
       iival = qrm_dunit
    else
       err = 23
       call qrm_error_print(err, name, aed=string)
       goto 9999
    end if

9999 continue
    if(present(info)) info = err
    return

  end subroutine qrm_ggetii




  !> @brief Used for counting the actual flops
  function qrm_count_realflops(m, n, k, op)
    real(kind(1.d0))                        :: qrm_count_realflops
    integer                                 :: m, k, n
    character                               :: op*(*)

    real(kind(1.d0))                        :: rk, rm, rn

    rk = real(k, kind(1.d0))
    rm = real(m, kind(1.d0))
    rn = real(n, kind(1.d0))

    qrm_count_realflops = 0.d0

    select case(op)
    case('panel')
       if(m .ge. k) then
          qrm_count_realflops = (2.d0*rk*rk*(rm-rk/3.d0))
       else
          rk = rm
          rn = real(n-m, kind(1.d0))
          qrm_count_realflops = (2.d0*rm*rm*(rm-rm/3.d0)) + (rk*rn*(4.d0*rm-rk))
       end if
    case('update')
       qrm_count_realflops = rk*rn*(4.d0*rm-rk)
    end select

    return
  end function qrm_count_realflops


  !> @brief Used for counting the real flops (i.e., it ignores the zeros
  !> included by the blocking)
  function qrm_count_pureflops(stair, n, j, nb)
    implicit none
    real(kind(1.d0))                        :: qrm_count_pureflops
    integer                                 :: stair(:)
    integer                                 :: n, nb, j

    integer                                 :: i

    qrm_count_pureflops=0
    do i=j, min(j+nb-1, n)
       qrm_count_pureflops = qrm_count_pureflops+(stair(i)-i+1)*(3 + 4*(n-i))
    end do

    return
  end function qrm_count_pureflops


end module qrm_common_mod
