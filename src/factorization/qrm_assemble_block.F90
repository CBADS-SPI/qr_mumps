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


subroutine _qrm_assemble_block(qrm_mat, fnum, br, bc)
  ! this subroutine performs the assembly of block (br,bc)
  ! of front fnum into its father
  !
  use _qrm_spmat_mod
  use _qrm_fdata_mod
  use qrm_common_mod
  use qrm_error_mod
  implicit none

  type(_qrm_spmat_type), target  :: qrm_mat
  integer                        :: fnum, br, bc

  type(_qrm_front_type), pointer :: front, father
  integer                        :: f, i, j, frow, fcol, fbc, fbr, fr, fc, m, n, npiv
  _qrm_data, pointer             :: blk(:,:)

  front  => qrm_mat%fdata%front_list(fnum)
  f      =  qrm_mat%adata%parent(fnum)
  father => qrm_mat%fdata%front_list(f)
  blk    => front%bc(br,bc)%c
  
  fr = (br-1)*front%mb - front%npiv
  fc = (bc-1)*front%nb - front%npiv

  m = size(blk,1)
  n = size(blk,2)

  do j=1, n
     if(fc+j .lt. 1) cycle ! still on the left of the CB
     fcol = front%colmap(fc+j)
     fbc  = (fcol-1)/father%nb+1
     fcol = mod(fcol-1,father%nb)+1
     do i=1, m
        if(fr+i .lt. 1) cycle ! still above the CB
        if(fr+i .gt. fc+j) exit ! below the diagonal, no need to keep going down
        frow = front%rowmap(fr+i)
        fbr  = (frow-1)/father%mb+1
        frow = mod(frow-1,father%mb)+1
        father%bc(fbr,fbc)%c(frow,fcol) = blk(i,j)
     end do
  end do

  return

end subroutine _qrm_assemble_block
