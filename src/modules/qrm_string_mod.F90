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
!> @file qrm_string_mod.F90
!! This file contains a module that implements string handling tools
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


!> @brief This module contains various string handling routines
module qrm_string_mod


  character(len=*), parameter   :: lcase='abcdefghijklmnopqrstuvwxyz'
  character(len=*), parameter   :: ucase='ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  interface qrm_str_tolower
     module procedure qrm_str_tolower
  end interface

  interface qrm_str_toupper
     module procedure qrm_str_toupper
  end interface qrm_str_toupper

  interface qrm_f2c_str
     module procedure qrm_f2c_str
  end interface qrm_f2c_str

contains 

  function  qrm_str_tolower(string)
    character(len=*), intent(in)  :: string
    character(len=len(string))    :: qrm_str_tolower
    integer  :: i,k
    
    do i=1,len(string)
       k = index(ucase, string(i:i))
       if (k /= 0) then 
          qrm_str_tolower(i:i) = lcase(k:k)
       else          
          qrm_str_tolower(i:i) = string(i:i)
       end if
    enddo
   
    return
    
  end function qrm_str_tolower
 
  function  qrm_str_toupper(string)
    character(len=*), intent(in)  :: string
    character(len=len(string))    :: qrm_str_toupper
    integer  :: i,k
    
    do i=1,len(string)
       k = index(lcase, string(i:i))
       if (k /= 0) then 
          qrm_str_toupper(i:i) = ucase(k:k)
       else          
          qrm_str_toupper(i:i) = string(i:i)
       end if
    enddo
   
    return
    
  end function qrm_str_toupper

  function qrm_f2c_str(strin) result(strout)
    use iso_c_binding
    character(len=*) :: strin
    character        :: strout(len(strin)+1)
    integer          :: i

    i = len(strin)
    strout(1:i+1) = transfer(strin(1:i)//c_null_char, strout(1:i+1))
    
    return
  end function qrm_f2c_str

  
end module qrm_string_mod
