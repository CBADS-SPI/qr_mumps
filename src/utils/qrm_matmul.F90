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
!> @file qrm_matmul.F90
!! This file contains a routine that does the sparse matrix - dense matrix  product 
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


!> @brief This subroutine computes the product y=beta*y+alpha*op(A)*x where
!! op(A) is either A or A' depending on the value of transp. 
!!
!! @param[in] qrm_mat the inpur A matrix
!!
!! @param[in] transp  if transp=_qrm_transp, op(A)=A'.
!!                    A otherwise
!!                     
!! @param[in] x       the source vector
!!                     
!! @param[in] alpha   a scalar
!!                     
!! @param[in] beta    a scalar
!!                     
!! @param[in,out] y   the destination vector          
!!
subroutine _qrm_matmul2d(qrm_mat, transp, alpha, x, beta, y)

  use _qrm_spmat_mod
  use qrm_string_mod
  implicit none

  type(_qrm_spmat_type)   :: qrm_mat
  _qrm_data, intent(out)  :: y(:,:)
  _qrm_data, intent(in)   :: x(:,:)
  _qrm_data, intent(in)   :: alpha, beta
  character(len=*)        :: transp

  integer                 :: nb, nrhs, j, c, r, k, i

  call qrm_get(qrm_mat, 'qrm_rhsnb', nb)
  nrhs = size(x,2)
  if(nb.le.0) nb = nrhs

  ! TODO: add checks for sizes etc.

  if(beta .eq. _qrm_zero) then
     y = _qrm_zero
  else
     y = beta*y
  end if

  ! shortcut
  if(alpha .eq. _qrm_zero) then
     return
  end if

  do k=1, nrhs, nb
     do i=1, qrm_mat%nz
        if((qrm_str_tolower(transp(1:1)) .eq. _qrm_transp) .or. (qrm_str_tolower(transp(1:1)) .eq. 't')) then
           c = qrm_mat%irn(i)
           r = qrm_mat%jcn(i)
           y(r,k:min(nrhs,k+nb-1)) = y(r,k:min(nrhs,k+nb-1))+alpha*_conjg(qrm_mat%val(i))*x(c,k:min(nrhs,k+nb-1))
        else
           r = qrm_mat%irn(i)
           c = qrm_mat%jcn(i)
           y(r,k:min(nrhs,k+nb-1)) = y(r,k:min(nrhs,k+nb-1))+alpha*qrm_mat%val(i)*x(c,k:min(nrhs,k+nb-1))
        end if
     end do
  end do

  return

end subroutine _qrm_matmul2d



!> @brief This subroutine computes the product y=beta*y+alpha*op(A)*x where
!! op(A) is either A or A' depending on the value of transp. This is
!! the vector version
!!
!! @param[in] qrm_mat the inpur A matrix
!!
!! @param[in] transp  if transp=_qrm_transp, op(A)=A'.
!!                    A otherwise
!!                     
!! @param[in] x       the source vector
!!                     
!! @param[in] alpha   a scalar
!!                     
!! @param[in] beta    a scalar
!!                     
!! @param[in,out] y   the destination vector          
!!
subroutine _qrm_matmul1d(qrm_mat, transp, alpha, x, beta, y)

  use _qrm_spmat_mod
  use qrm_string_mod
  implicit none

  type(_qrm_spmat_type)   :: qrm_mat
  _qrm_data, intent(out)  :: y(:)
  _qrm_data, intent(in)   :: x(:)
  _qrm_data, intent(in)   :: alpha, beta
  character(len=*)        :: transp

  integer                 :: c, r, i, n

  n = size(x,1)

  ! TODO: add checks for sizes etc.

  if(beta .eq. _qrm_zero) then
     y = _qrm_zero
  else
     y = beta*y
  end if

  ! shortcut
  if(alpha .eq. _qrm_zero) then
     return
  end if

  do i=1, qrm_mat%nz
     if((qrm_str_tolower(transp(1:1)) .eq. _qrm_transp) .or. (qrm_str_tolower(transp(1:1)) .eq. 'c')) then
        c = qrm_mat%irn(i)
        r = qrm_mat%jcn(i)
        y(r) = y(r)+alpha*_conjg(qrm_mat%val(i))*x(c)
     else
        r = qrm_mat%irn(i)
        c = qrm_mat%jcn(i)
        y(r) = y(r)+alpha*qrm_mat%val(i)*x(c)
     end if
  end do

  return

end subroutine _qrm_matmul1d
