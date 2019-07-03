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
!> @file qrm_mod.F90
!! This file contains the module which groups all the qr_mumps useful modules.
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


!> This is intended to be the only module to be included in user's
!! programs
module _qrm_mod

  use _qrm_spmat_mod
  use _qrm_utils_mod
  use qrm_mem_mod
  use qrm_error_mod
  use _qrm_analysis_mod
  use _qrm_factorization_mod
  use _qrm_solve_mod
  use _qrm_sdata_mod
  use qrm_common_mod
  use _qrm_methods_mod
  use qrm_const_mod
  use qrm_dscr_mod
  use qrm_string_mod

end module _qrm_mod



