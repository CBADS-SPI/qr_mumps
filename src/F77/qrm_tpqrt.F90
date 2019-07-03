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

! This routine is a variant of the LAPACK _TPQRT routine which allows
! for handling matrices that have a staircase structure

subroutine _qrm_tpqrt( m, n, l, nb, stair, ofsa, ofsb, a, lda, b, ldb, t, ldt, work, info )

  integer   :: info, lda, ldb, ldt, n, m, l, nb, ofsa, ofsb
  integer   :: stair(*)
  _qrm_data :: a( lda, * ), b( ldb, * ), t( ldt, * ), work( * )

  integer   ::  i, j, ib, im, in, il

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !   the old working version starts here
  !      i = 1
  !      do j=1, n
  !         if(stair(j) .gt. ofsb) exit
  !         if(stair(j) .gt. ofsa) i = i+1
  !      end do
  !      if(j.gt.n) return

  !      in = n-j+1
  !      if(l.gt.0) then
  !         il = min(in, stair(j+in-1)-ofsb)
  !         im = il
  !      else
  !         im = min(m,stair(j+in-1)-ofsb)
  !         il = l
  !      end if
  !      ib = min(in,nb)

  !c      write(*,*)'tpqrt  ',im,in,il
  !      call dtpqrt( im, in, il, ib, a(i,j), lda,
  !     $     b(1,j), ldb, t(1,j), ldt, work, info )
  !   the old working version ends here
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


  if((l.ne.0) .and. (l.ne.m)) then
     call xerbla( '_xtpqrt', -9999 )
     return
  end if

  i = 1
  do js=1, n
     if(stair(js) .gt. ofsb) exit
     if(stair(js) .gt. ofsa) i = i+1
  end do

  if(js.gt.n) return

  do j=js, n, nb
     ib = min( n-j+1, nb )
     if(l.eq.0) then
        ! todo: check
        im = min(stair(j+ib-1)-ofsb,m)
        il = 0
     else if (l.eq.m) then
        im = min(j+ib-js,m)
        !            il = min(min(im,ib),m-j+js)
        if(j-js.gt.l) then
           il = 0
        else
           il = im - j + js
        end if
     end if

     call _xtpqrt2( im, ib, il, a(i,j), lda, b( 1, j ), ldb, &
          & t(1, j ), ldt, iinfo )

     if( j+ib.le.n ) then
        call _xtprfb( 'l', _qrm_transp, 'f', 'c', im, n-j-ib+1, &
             & ib, il, b( 1, j ), ldb, t( 1, j ), ldt, &
             & a( i, j+ib ), lda, b( 1, j+ib ), ldb, &
             & work, ib )
     end if

     i = i+ib
  end do

  return

end subroutine _qrm_tpqrt
