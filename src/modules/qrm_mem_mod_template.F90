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


!beg_template: alloc
subroutine qrm_#pref#alloc_#rank##arith#(a, #sizes#, info)

  #type##prec#, #array#, dimension(#dim#) :: a
  integer, intent(in)                     :: #sizes#
  integer, optional                       :: info

  integer :: err

  err = 0
  
  if(#checksize#) return

  if(#check#(a)) then
     err = 4
  else
     allocate(a(#sizes#), stat=err)
     if(err .eq. 0) then
        call qrm_mem_upd(size(a,kind=qrm_mem_kind)*qrm_sizeof_#arith#_)
     else
        err= 12
     end if
  end if

  if(present(info)) info = err
  return

end subroutine qrm_#pref#alloc_#rank##arith#
!end_template: alloc


!beg_template: dealloc
subroutine qrm_#pref#dealloc_#rank##arith#(a, info)

  #type##prec#, #array#, dimension(#dim#) :: a
  integer, optional                       :: info

  integer(kind=qrm_mem_kind) :: s
  integer :: err

  err = 0
  
  if(#check#(a)) then
     s = size(a,kind=qrm_mem_kind)
     deallocate(a, stat=err)
  else
     return
  end if
  if(err .eq. 0) then
     call qrm_mem_upd(-s*qrm_sizeof_#arith#_)
  else
     err = 6
  end if

  if(present(info)) info = err
  return

end subroutine qrm_#pref#dealloc_#rank##arith#
!end_template: dealloc

!beg_template: size
function qrm_#pref#size_#rank##arith#(a)

  implicit none

  integer :: qrm_#pref#size_#rank##arith#
  #type##prec#, #array#, dimension(#dim#) :: a

  if(#check#(a)) then
     ! a is allocated    
     qrm_#pref#size_#rank##arith# = size(a)
  else
     qrm_#pref#size_#rank##arith# = 0
  end if

  return

end function qrm_#pref#size_#rank##arith#
!end_template: size


!beg_template: move_alloc
subroutine qrm_#pref#move_alloc_#rank##arith#(src, dst)
  #type##prec#, #array#, dimension(#dim#) :: src, dst
  #move#
  return
end subroutine qrm_#pref#move_alloc_#rank##arith#
!end_template: move_alloc


!beg_template: realloc
subroutine qrm_#pref#realloc_#rank##arith#(a, n, info, copy)

  #type##prec#, #array#, dimension(#dim#) :: a
  integer                                 :: n
  logical, optional                       :: copy
  integer, optional                       :: info

  integer                                 :: i, err
  logical                                 :: icopy
  #type##prec#, #array#, dimension(#dim#) :: tmp

  ! a is not allocated. just allocate it and return
  if(.not.#check#(a)) then
     call qrm_#pref#alloc_#rank##arith#(a, n, err)
     goto 9999
  else
     if(size(a,kind=qrm_mem_kind).ge.n) return
  end if

  if(present(copy)) then
     icopy=copy
  else
     icopy = .false.
  end if

  ! we need to reallocate
  if(icopy) then
     ! we must save a copy
     call qrm_#pref#move_alloc_#rank##arith#(a, tmp)
  else
     call qrm_#pref#dealloc_#rank##arith#(a)
  end if

  call qrm_#pref#alloc_#rank##arith#(a, n, err)
  if(err.ne.0) goto 9999

  ! check if copy is to be done
  if(icopy) then
     do i=1, size(a)
        a(i) = tmp(i)
     end do
     call qrm_#pref#dealloc_#rank##arith#(tmp, err)
  end if

9999 continue
  if(present(info)) info = err
  return
  
  return

end subroutine qrm_#pref#realloc_#rank##arith#
!end_template: realloc


