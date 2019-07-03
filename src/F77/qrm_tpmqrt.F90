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

! This routine is a variant of the LAPACK _TPMQRT routine which allows
! for handling matrices that have a staircase structure

subroutine _qrm_tpmqrt( side, trans, m, n, k, l, nb, stair, &
     & ofsa, ofsb, v, ldv, t, ldt, &
     & a, lda, b, ldb, work, info )

  character :: side, trans
  integer   :: info, k, ldv, lda, ldb, m, n, l, nb, ldt, ofsa, ofsb
  integer   :: stair(*)
  _qrm_data :: v( ldv, * ), a( lda, * ), b( ldb, * ), &
       & t( ldt, * ), work( * )

  integer   ::  i, j, ib, im, in, il, ik 
  logical   ::  lsame

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !   the old working version starts here
  !      i = 1
  !      do j=1, k
  !         if(stair(j) .gt. ofsb) exit
  !         if(stair(j) .gt. ofsa) i = i+1
  !      end do
  !      if(j.gt.k) return
  !
  !      
  !      ik = k-j+1
  !      if(l.gt.0) then
  !         il = min(ik, stair(j+ik-1)-ofsb)
  !         im = il
  !      else
  !         im = min(m,stair(j+ik-1)-ofsb)
  !         il = l
  !      end if
  !      ib = min(ik,nb)
  !
  !      call dtpmqrt( side, trans, im, n, ik, il, ib,
  !     $     v(1,j), ldv, t(1,j), ldt,
  !     $     a(i,1), lda, b(1,1), ldb, work, info )
  !   the old working version ends here
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  if((l.ne.0) .and. (l.ne.m)) then
     call xerbla( '_xtpmqrt', -9999 )
     return
  end if

  if(.not. lsame(side,'l')) then
     call xerbla( '_xtpmqrt', -9998 )
     return
  end if

  is = 1
  do js=1, k
     if(stair(js) .gt. ofsb) exit
     if(stair(js) .gt. ofsa) is = is+1
  end do

  if(js.gt.k) return

  if(lsame(trans,_qrm_transp)) then
     i = is
     do j=js, k, nb
        ib = min( k-j+1, nb )
        if(l.eq.0) then
           ! todo: check
           im = min(stair(j+ib-1)-ofsb,m)
           il = 0
        else if (l.eq.m) then
           im = min(j+ib-js,m)
           if(j-js.gt.l) then
              il = 0
           else
              il = im - j + js
           end if
        end if
        call _xtprfb( 'l', _qrm_transp, 'f', 'c', im, n, ib, il, &
             & v( 1, j ), ldv, t( 1, j ), ldt, &
             & a( i, 1 ), lda, b, ldb, work, ib )

        i = i+ib
     end do
  else if(lsame(trans,'n')) then

     je = js + ((k-js)/nb)*nb
     i  = is + ((k-js)/nb)*nb
     do j=je, js, -nb
        ib = min( k-j+1, nb )
        if(l.eq.0) then
           ! todo: check
           im = min(stair(j+ib-1)-ofsb,m)
           il = 0
        else if (l.eq.m) then
           im = min(j+ib-js,m)
           if(j-js.gt.l) then
              il = 0
           else
              il = im - j + js
           end if
        end if
        call _xtprfb( 'l', 'n', 'f', 'c', im, n, ib, il, &
             & v( 1, j ), ldv, t( 1, j ), ldt, &
             & a( i, 1 ), lda, b, ldb, work, ib )

        i = i-nb
     end do


  end if

  return

end subroutine _qrm_tpmqrt
