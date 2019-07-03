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
!> @file qrm_solve_sing_front.F90
!! This file contains a routine that handles the singletons front during the solve
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief This function handles the front containing the singletons
!! during the solve for R or R'
!! 
!! @param[in] qrm_mat   the main qrm data structure after factorization.
!!
!! @param[in]     b     a 2d array containing the RHS vectors
!!
!! @param[out]    x     a 2d array containing the solution vectors
!!
!! @param[in] trans     a string saying whether R or R^T will be solved
!!                      for. Only the first character is important.
!! 
subroutine _qrm_solve_sing_front(qrm_mat, b, x, trans)

  use _qrm_spmat_mod
  use qrm_string_mod
  use _qrm_fdata_mod
  implicit none

  type(_qrm_spmat_type), target :: qrm_mat
  _qrm_data, intent(inout)      :: b(:,:)
  _qrm_data, intent(inout)      :: x(:,:)
  character                      :: trans

  type(_qrm_front_type), pointer :: front
  integer    :: i, j, row, col
  _qrm_data :: d

  ! the front with the singletons is always the first one
  front => qrm_mat%fdata%front_list(1)

  ! do the computations.
  !
  ! THIS ROUTINE ASSUMES THAT THE DIAGONAL ELEMENTS IN FRONT%AVAL AND
  ! FRONT%AJCN ALWAYS COME BEFORE THE OTHERS ALONG THEIR ROW. THIS IS
  ! ENFORCED INSIDE QRM_FACTORIZATION_INIT
  !
  if(qrm_str_tolower(trans) .eq. 'n') then
     do i=front%anrows, 1, -1
        row =qrm_mat%adata%rperm(i)
        do j=front%aiptr(i+1)-1, front%aiptr(i)+1,-1
           col = front%ajcn(j)
           b(row,:) = b(row,:) - front%aval(j)*x(col,:)
        end do
        x(qrm_mat%adata%cperm(i),:) = b(row,:)/front%aval(front%aiptr(i))
     end do

  else if(qrm_str_tolower(trans) .eq. _qrm_transp) then

     do j=1, front%anrows
        col = qrm_mat%adata%rperm(j)
        x(col,:) = b(qrm_mat%adata%cperm(j),:)/_conjg(front%aval(front%aiptr(j)))
        do i=front%aiptr(j)+1, front%aiptr(j+1)-1
           row = front%ajcn(i)
           b(row,:) = b(row,:)-_conjg(front%aval(i))*x(col,:)
        end do
     end do
  else

  end if

  return

end subroutine _qrm_solve_sing_front

