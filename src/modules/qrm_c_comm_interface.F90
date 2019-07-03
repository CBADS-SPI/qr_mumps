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
!> @file qrm_cintface.F90
!! This file contains the C interface for qr_mumps.
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief This module contains the definition of the qr_mumps C interface common
!! to all precisions/types.
module qrm_c_comm_interface
  use iso_c_binding
  use qrm_error_mod
  use qrm_common_mod
  
contains

  !> @brief C equivalent of the @link _qrm_dscr_mod::qrm_init @endlink
  !> routine
  function qrm_init_c(nthreads) result(info) bind(c)
    use qrm_dscr_mod
    implicit none
    integer(c_int), value  :: nthreads
    integer(c_int)         :: info

    if(nthreads.gt.0) then
       call qrm_init(nthreads, info)
    else
       call qrm_init(info=info)
    end if

    return
  end function qrm_init_c
    
  !> @brief C equivalent of the @link _qrm_dscr_mod::qrm_finalize @endlink
  !> routine
  subroutine qrm_finalize_c() bind(c)
    use qrm_dscr_mod
    implicit none

    call qrm_finalize()

    return
  end subroutine qrm_finalize_c
    
  !> @brief C equivalent of the @link _qrm_common_mod::_qrm_gseti @endlink
  !> routine (only for global, integer type)
  function qrm_gseti_c(string, val) result(info) bind(c)
    character(kind=c_char) :: string(40)
    integer(c_int), value  :: val
    integer(c_int)         :: info
    
    character(len=40) :: a
    
    write(a,'(40a)')string
    
    call qrm_gseti(a, val, info)

    return
    
  end function qrm_gseti_c

  !> @brief C equivalent of the @link _qrm_common_mod::_qrm_ggeti @endlink
  !> routine (only for global, integer type)
  function qrm_ggeti_c(string, val) result(info) bind(c)
    character(kind=c_char) :: string(40)
    integer(c_int)         :: val
    integer(c_int)         :: info
    
    character(len=40) :: a
    
    write(a,'(40a)')string
    
    call qrm_ggeti(a, val, info)

    return
    
  end function qrm_ggeti_c

  !> @brief C equivalent of the @link _qrm_common_mod::_qrm_ggetii @endlink
  !> routine (only for global, integer type)
  function qrm_ggetii_c(string, val) result(info) bind(c)
    character(kind=c_char) :: string(40)
    integer(c_long_long)   :: val
    integer(c_int)         :: info
    
    character(len=40) :: a
    
    write(a,'(40a)')string
    
    call qrm_ggetii(a, val, info)

    return
    
  end function qrm_ggetii_c

end module qrm_c_comm_interface

