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
!> @file qrm_error_mod.F90
!! This file contains the module that implements the error management
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief This module contains all the error management routines and 
!! data. 


module qrm_error_mod
  use qrm_const_mod

contains

    subroutine qrm_error_print(code, where, ied, aed)
    implicit none

    integer                         :: code
    character(len=*)                :: where
    integer, optional               :: ied(:)
    character(len=*), optional      :: aed

    
    __QRM_PRNT_ERR('("Error in subroutine ",a30, " :")')where
    select case(code)
    case(0)
       __QRM_PRNT_ERR('("Generic error")')
    case(1)
       __QRM_PRNT_ERR('("Sparse matrix format ",a3," is not (yet) supported.")')aed(1:3)
    case(2)
       __QRM_PRNT_ERR('("Symmetric matrices are not supported.")')
    case(3)
       __QRM_PRNT_ERR('("qrm_spmat%cntl is not associated/valid.")')
    case(4)
       __QRM_PRNT_ERR('("Trying to allocate an already allocated array.")')
    case(5)
       __QRM_PRNT_ERR('("Memory allocation problem. Size required: ",i30)')ied(1)
    case(6)
       __QRM_PRNT_ERR('("Memory deallocation problem.")')
    case(8)
       __QRM_PRNT_ERR('("Input column permutation not provided/valid")')
    case(9)
       __QRM_PRNT_ERR('("Requested ordering method unknown: ",i3)')ied(1)
    case(10)
       __QRM_PRNT_ERR('("Insufficient size for array: ",a20)')aed(1:20)
    case(11)
       __QRM_PRNT_ERR('("Error in lapack routine: ",i3)')ied(1)
    case(12)
       __QRM_PRNT_ERR('("Out of memory")')
    case(13)
       __QRM_PRNT_ERR('("The analysis must be done before the factorization")')
    case(14)
       __QRM_PRNT_ERR('("The factorization must be done before the solve")')
    case(15)
       __QRM_PRNT_ERR('("This type of norm is not implemented.")')
    case(16)
       __QRM_PRNT_ERR('("Requested ordering method not available: ",a20)')aed
    case(17)
       __QRM_PRNT_ERR('("Error from call to subroutine ",a30,": ",i3)')aed,ied(1)
    case(18)
       __QRM_PRNT_ERR('("COLAMD error ")')
    case(19)
       __QRM_PRNT_ERR('("SCOTCH error ")')
    case(20)
       __QRM_PRNT_ERR('("Factorization error ")')
    case(21)
       __QRM_PRNT_ERR('("Apply error ")')
    case(22)
       __QRM_PRNT_ERR('("Solve error ")')
    case(23)
       __QRM_PRNT_ERR('("Incorrect argument to qrm_set/qrm_get ",a30)')aed
    case(25)
       __QRM_PRNT_ERR('("Problem opening file ",a30)')aed
    case(26)
       __QRM_PRNT_ERR('("Unknown error action ",i10)')ied(1)
    case(27)
       __QRM_PRNT_ERR('("Incompatible values in qrm_spmat%icntl ",i2,2x,i2)')ied(1:2)
    case(28)
       __QRM_PRNT_ERR('("Incorrect value for qrm_mb_/qrm_nb_/qrm_ib_ :",i4,2x,i4,2x,i4)')ied(1:3)
    case(29)
       __QRM_PRNT_ERR('("Incorrect value for qrm_spmat%m/n/nz ",i6,2x,i6,2x,i16)')ied(1:3)
    case(30)
       __QRM_PRNT_ERR('("qrm_apply cannot be called if the H matrix is discarded.")')
    case(31)
       __QRM_PRNT_ERR('("StarPU initialization error.")')
    case default
       __QRM_PRNT_ERR('("Unknown error code",i4)')code
    end select

    __QRM_PRNT_ERR('(" ")')
    return

  end subroutine qrm_error_print


end module qrm_error_mod
