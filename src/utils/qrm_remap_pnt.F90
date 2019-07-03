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
!> @file qrm_remap_pnt.F90
!! This file contains a routine that does a dirty trick to convert a 1d array into a 2d array.
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


!> @brief This function makes a 2D pointer point to
!! a 1D array. 

!> This is needed to have a single code handling single
!! and multiple right-hand-sides. This is a dirty workaround that
!! relies on iso_c_binding while waiting for the array bounds
!! remapping to be supported by compilers (in F2008 standard).
!! 
!! *Input*:
!! arr1d - the input 1d array
!! 
!! *Output*:
!! pnt2d - the output 2d pointer pointing to arr1d
subroutine _qrm_remap_pnt(arr1d, pnt2d, n)

  use iso_c_binding
  implicit none

  integer :: n
  _qrm_data, target  :: arr1d(1:n)
  _qrm_data, pointer :: pnt2d(:,:)

  type(c_ptr) :: cptr

  cptr = c_loc(arr1d(1))

  call c_f_pointer(cptr, pnt2d, (/n,1/))

  return

end subroutine _qrm_remap_pnt

