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


!!##############################################################################################
!> @file qrm_prnt_array.F90
!! FIXME: add comments
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!!##############################################################################################


subroutine qrm_prnt_iarray(a, lab, unit)
  ! Function: qrm_prnt_iarray
  !
  ! *Input*:
  !
  ! *Output*:
  !

  integer   :: a(:)
  character :: lab*(*)
  integer, optional :: unit
  integer   :: mn, mx, s, i, lg, iunit
  character :: fmt*10, fmt2*10


  if(present(unit)) then
     iunit = unit
  else
     iunit = 6
  end if

  write(iunit,'(a15,"= [ ")',advance='no')lab
  do i=1, size(a)
     mx = abs(a(i))
     lg = 10
     s  = 1
     do
        if (lg .gt. mx) exit
        s = s+1
        lg=lg*10
     end do
     if(a(i) .lt. 0) s = s+1
     fmt=''
     write(fmt,'(i4)')s
     fmt=adjustl(fmt)
     write(fmt2,'("(i",a4,",1x)")')fmt(1:4)
     write(iunit,fmt2,advance='no')a(i)
  end do
  
  write(iunit,'(" ];")')
  write(iunit,'(" ")')
   
  return

end subroutine qrm_prnt_iarray


subroutine qrm_prnt_sarray(a, lab, unit)
  ! Function: qrm_prnt_iarray
  !
  ! *Input*:
  !
  ! *Output*:
  !

  real(kind(1.e0))  :: a(:)
  character         :: lab*(*)
  integer, optional :: unit

  integer   :: mn, mx, s, i, lg, iunit
  character :: fmt*12, fmt2*12


  if(present(unit)) then
     iunit = unit
  else
     iunit = 6
  end if

  write(iunit,'(a15,"= [ ")',advance='no')lab
  do i=1, size(a)
     mx = floor(abs(a(i)))
     lg = 10
     s  = 1
     do
        if (lg .gt. mx) exit
        s = s+1
        lg=lg*10
     end do
     if(a(i) .lt. 0) s = s+1
     fmt=''
     write(fmt,'(i4,".4")')s+5
     fmt=adjustl(fmt)
     write(fmt2,'("(f",a5,",1x)")')fmt(1:5)
     write(iunit,fmt2,advance='no')a(i)
  end do
  
  write(iunit,'(" ];")')
  write(iunit,'(" ")')
   
  return

end subroutine qrm_prnt_sarray



subroutine qrm_prnt_darray(a, lab, unit)
  ! Function: qrm_prnt_iarray
  !
  ! *Input*:
  !
  ! *Output*:
  !

  real(kind(1.d0))  :: a(:)
  character         :: lab*(*)
  integer, optional :: unit

  integer   :: mn, mx, s, i, lg, iunit
  character :: fmt*12, fmt2*12


  if(present(unit)) then
     iunit = unit
  else
     iunit = 6
  end if

  write(iunit,'(a15,"= [ ")',advance='no')lab
  do i=1, size(a)
     mx = floor(abs(a(i)))
     lg = 10
     s  = 1
     do
        if (lg .gt. mx) exit
        s = s+1
        lg=lg*10
     end do
     if(a(i) .lt. 0) s = s+1
     fmt=''
     write(fmt,'(i4,".14")')s+15
     fmt=adjustl(fmt)
     write(fmt2,'("(f",a6,",1x)")')fmt(1:6)
     write(iunit,fmt2,advance='no')a(i)
  end do
  
  write(iunit,'(" ];")')
  write(iunit,'(" ")')
   
  return

end subroutine qrm_prnt_darray


subroutine qrm_prnt_d2array(a, lab, unit)
  ! Function: qrm_prnt_iarray
  !
  ! *Input*:
  !
  ! *Output*:
  !
  implicit none
  
  real(kind(1.d0))  :: a(:,:)
  character         :: lab*(*)
  integer, optional :: unit

  integer   :: mn, mx, s, i, j, lg, iunit
  character :: fmt*12, fmt2*12


  if(present(unit)) then
     iunit = unit
  else
     iunit = 6
  end if

  write(iunit,'(a15,"= [ ")',advance='no')lab
  do i=1, size(a,1)
     do j=1, size(a,2)
        mx = floor(abs(a(i,j)))
        lg = 10
        s  = 1
        do
           if (lg .gt. mx) exit
           s = s+1
           lg=lg*10
        end do
        if(a(i,j) .lt. 0) s = s+1
        fmt=''
        write(fmt,'(i4,".14")')s+15
        fmt=adjustl(fmt)
        write(fmt2,'("(f",a6,",1x)")')fmt(1:6)
        write(iunit,fmt2,advance='no')a(i,j)
     end do
     write(iunit,'(" ")')
  end do
  write(iunit,'(" ];")')
  write(iunit,'(" ")')
   
  return
  
end subroutine qrm_prnt_d2array






subroutine qrm_prnt_carray(a, lab, unit)
  ! Function: qrm_prnt_carray
  !
  ! *Input*:
  !
  ! *Output*:
  !

  complex(kind(1.e0))  :: a(:)
  character         :: lab*(*)
  integer, optional :: unit

  integer   :: mn, mx, s, i, lg, iunit
  character :: fmt*12, fmt2*12, fmt3*12


  if(present(unit)) then
     iunit = unit
  else
     iunit = 6
  end if

  write(iunit,'(a15,"= [ ")',advance='no')lab
  do i=1, size(a)
     mx = floor(abs(real(a(i))))
     lg = 10
     s  = 1
     do
        if (lg .gt. mx) exit
        s = s+1
        lg=lg*10
     end do
     if(real(a(i)) .lt. 0) s = s+1
     fmt=''
     write(fmt,'(i4,".4")')s+5
     fmt=adjustl(fmt)
     write(fmt2,'("(f",a5,","","",")')fmt(1:5)

     mx = floor(abs(aimag(a(i))))
     lg = 10
     s  = 1
     do
        if (lg .gt. mx) exit
        s = s+1
        lg=lg*10
     end do
     if(aimag(a(i)) .lt. 0) s = s+1
     fmt=''
     write(fmt,'(i4,".4")')s+5
     fmt=adjustl(fmt)
     write(fmt3,'("f",a5,",1x)")')fmt(1:5)

     write(iunit,fmt2//fmt3,advance='no')a(i)
  end do
  
  write(iunit,'(" ];")')
  write(iunit,'(" ")')
   
  return

end subroutine qrm_prnt_carray



subroutine qrm_prnt_zarray(a, lab, unit)
  ! Function: qrm_prnt_zarray
  !
  ! *Input*:
  !
  ! *Output*:
  !

  complex(kind(1.d0))  :: a(:)
  character         :: lab*(*)
  integer, optional :: unit

  integer   :: mn, mx, s, i, lg, iunit
  character :: fmt*12, fmt2*12, fmt3*12


  if(present(unit)) then
     iunit = unit
  else
     iunit = 6
  end if

  write(iunit,'(a15,"= [ ")',advance='no')lab
  do i=1, size(a)
     mx = floor(abs(real(a(i))))
     lg = 10
     s  = 1
     do
        if (lg .gt. mx) exit
        s = s+1
        lg=lg*10
     end do
     if(real(a(i)) .lt. 0) s = s+1
     fmt=''
     write(fmt,'(i4,".4")')s+5
     fmt=adjustl(fmt)
     write(fmt2,'("(f",a5,","","",")')fmt(1:5)

     mx = floor(abs(aimag(a(i))))
     lg = 10
     s  = 1
     do
        if (lg .gt. mx) exit
        s = s+1
        lg=lg*10
     end do
     if(aimag(a(i)) .lt. 0) s = s+1
     fmt=''
     write(fmt,'(i4,".4")')s+5
     fmt=adjustl(fmt)
     write(fmt3,'("f",a5,",1x)")')fmt(1:5)

     write(iunit,fmt2//fmt3,advance='no')a(i)
  end do
  
  write(iunit,'(" ];")')
  write(iunit,'(" ")')
   
  return

end subroutine qrm_prnt_zarray
