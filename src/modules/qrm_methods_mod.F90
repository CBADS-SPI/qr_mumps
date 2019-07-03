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
!> @file qrm_methods_mod.F90
!! This file contains a module with all the generic interfaces for the typed methods routines.
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


!> @brief This module contains generic methods
module _qrm_methods_mod
  

  interface qrm_least_squares
     procedure               ::  _qrm_least_squares2d, _qrm_least_squares1d
  end interface qrm_least_squares

  interface _qrm_least_squares
     subroutine _qrm_least_squares2d(qrm_mat, b, x, info)
       use _qrm_spmat_mod
       type(_qrm_spmat_type) :: qrm_mat
       _qrm_data             :: b(:,:), x(:,:)
       integer, optional     :: info
     end subroutine _qrm_least_squares2d
     subroutine _qrm_least_squares1d(qrm_mat, b, x, info)
       use _qrm_spmat_mod
       type(_qrm_spmat_type) :: qrm_mat
       _qrm_data             :: b(:), x(:)
       integer, optional     :: info
     end subroutine _qrm_least_squares1d
  end interface _qrm_least_squares
  


  interface qrm_min_norm
     procedure               :: _qrm_min_norm2d, _qrm_min_norm1d
  end interface qrm_min_norm

  interface _qrm_min_norm
     subroutine _qrm_min_norm2d(qrm_mat, b, x, info)
       use _qrm_spmat_mod
       type(_qrm_spmat_type) :: qrm_mat
       _qrm_data             :: b(:,:), x(:,:)
       integer, optional     :: info
     end subroutine _qrm_min_norm2d
     subroutine _qrm_min_norm1d(qrm_mat, b, x, info)
       use _qrm_spmat_mod
       type(_qrm_spmat_type) :: qrm_mat
       _qrm_data             :: b(:), x(:)
       integer, optional     :: info
     end subroutine _qrm_min_norm1d
  end interface _qrm_min_norm
  

  interface qrm_residual_norm
     procedure               :: _qrm_residual_norm2d, _qrm_residual_norm1d
  end interface qrm_residual_norm

  interface _qrm_residual_norm
     subroutine _qrm_residual_norm2d(qrm_mat, b, x, nrm, info)
       use _qrm_spmat_mod
       _qrm_real             :: nrm(:)
       type(_qrm_spmat_type) :: qrm_mat
       _qrm_data             :: b(:,:), x(:,:)
       integer, optional     :: info
     end subroutine _qrm_residual_norm2d
     subroutine _qrm_residual_norm1d(qrm_mat, b, x, nrm, info)
       use _qrm_spmat_mod
       _qrm_real             :: nrm
       type(_qrm_spmat_type) :: qrm_mat
       _qrm_data             :: b(:), x(:)
       integer, optional     :: info
     end subroutine _qrm_residual_norm1d
  end interface _qrm_residual_norm



  interface qrm_residual_orth
     procedure               :: _qrm_residual_orth2d, _qrm_residual_orth1d
     procedure               :: _qrm_residual_and_orth2d, _qrm_residual_and_orth1d
  end interface qrm_residual_orth

  interface _qrm_residual_orth
     subroutine _qrm_residual_orth2d(qrm_mat, r, nrm, info)
       use _qrm_spmat_mod
       _qrm_real             :: nrm(:)
       type(_qrm_spmat_type) :: qrm_mat
       _qrm_data             :: r(:,:)
       integer, optional     :: info
     end subroutine _qrm_residual_orth2d
     subroutine _qrm_residual_orth1d(qrm_mat, r, nrm, info)
       use _qrm_spmat_mod
       _qrm_real             :: nrm
       type(_qrm_spmat_type) :: qrm_mat
       _qrm_data             :: r(:)
       integer, optional     :: info
     end subroutine _qrm_residual_orth1d
     subroutine _qrm_residual_and_orth2d(qrm_mat, b, x, nrm, info)
       use _qrm_spmat_mod
       _qrm_real             :: nrm(:)
       type(_qrm_spmat_type) :: qrm_mat
       _qrm_data             :: b(:,:), x(:,:)
       integer, optional     :: info
     end subroutine _qrm_residual_and_orth2d
     subroutine _qrm_residual_and_orth1d(qrm_mat, b, x, nrm, info)
       use _qrm_spmat_mod
       _qrm_real             :: nrm
       type(_qrm_spmat_type) :: qrm_mat
       _qrm_data             :: b(:), x(:)
       integer, optional     :: info
     end subroutine _qrm_residual_and_orth1d
  end interface _qrm_residual_orth


  interface qrm_get_r
     subroutine _qrm_get_r(qrm_mat, r, info)
       use _qrm_spmat_mod
       type(_qrm_spmat_type), target  :: qrm_mat
       type(_qrm_spmat_type)          :: r
       integer, optional              :: info
     end subroutine _qrm_get_r
  end interface qrm_get_r
  
end module _qrm_methods_mod
