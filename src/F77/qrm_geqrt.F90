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

! This routine is a variant of the LAPACK _GEQRT routine which allows
! for handling matrices that have a staircase structure

subroutine _qrm_geqrt( m, n, nb, stair, ofs, a, lda, t, ldt, work, info )

  implicit none

  integer            :: info, lda, ldt, m, n, nb, ofs

  integer            :: stair(n)
  _qrm_data          :: a( lda, * ), t( ldt, * ), work( * )

  integer            ::   i, j, ib, iinfo, k, im, st, en
  logical, parameter ::   use_recursive_qr=.true.

  info = 0
  if( m.lt.0 ) then
     info = -1
  else if( n.lt.0 ) then
     info = -2
  else if( nb.lt.1 .or. ( nb.gt.min(m,n) .and. min(m,n).gt.0 ) )then
     info = -3
  else if( lda.lt.max( 1, m ) ) then
     info = -5
  else if( ldt.lt.nb ) then
     info = -7
  end if
  if( info.ne.0 ) then
     call xerbla( '_geqrt', -info )
     return
  end if

  k = min( m, n )
  if( k.eq.0 ) return

  st = 1
  do while(stair(st)-ofs.le.0)
     st=st+1
     if(st.gt.n) return
  end do
  en = min(st+min(m,stair(n)-ofs)-1,n)

  i = 1
  do j = st, en,  nb
     ib = min( en-j+1, nb )
     im = max(min(stair(j+ib-1)-ofs,m) - i + 1,ib)
     if(im.gt.0) then

        if( use_recursive_qr ) then
           call _xgeqrt3( im, ib, a(i,j), lda, t(1,j), ldt, iinfo )
        else
           call _xgeqrt2( im, ib, a(i,j), lda, t(1,j), ldt, iinfo )
        end if

        if( j+ib.le.n ) then

           call _xlarfb('l', _qrm_transp, 'f', 'c', im, n-j-ib+1, ib, a(i,j), &
                     & lda, t(1,j), ldt, a(i,j+ib), lda, work , n )
        end if
     end if
     i = i+ib
  end do
  return
end subroutine _qrm_geqrt
