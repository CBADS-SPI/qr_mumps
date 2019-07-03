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


! This routine is a variant of the LAPACK _GEMQRT routine which allows
! for handling matrices that have a staircase structure

subroutine _qrm_gemqrt(side, trans,   &
     & m, n, k, nb, stair, ofs,       &
     & v, ldv,                        &
     & t, ldt,                        &
     & c, ldc,                        &
     & work, info)

  implicit none

  character :: side, trans
  integer   :: info, k, ldv, ldc, m, n, nb, ldt, ofs
  integer   :: stair(*)
  _qrm_data ::  v( ldv, * ), c( ldc, * ), t( ldt, * ), work( * )

  logical   ::         left, right, tran, notran
  integer   ::         i, j, ib, ldwork, kf, q, st, en, im
  logical   ::         lsame

  info   = 0
  left   = lsame( side,  'l' )
  right  = lsame( side,  'r' )
  tran   = lsame( trans, _qrm_transp )
  notran = lsame( trans, 'n' )

  if( left ) then
     ldwork = max( 1, n )
     q = m
  else if ( right ) then
     ldwork = max( 1, m )
     q = n
  end if
  if( .not.left .and. .not.right ) then
     info = -1
  else if( .not.tran .and. .not.notran ) then
     info = -2
  else if( m.lt.0 ) then
     info = -3
  else if( n.lt.0 ) then
     info = -4
     !      else if( k.lt.0 .or. k.gt.q ) then
     !         info = -5
  else if( nb.lt.1 .or. (nb.gt.k .and. k.gt.0)) then
     info = -6
  else if( ldv.lt.max( 1, q ) ) then
     info = -8
  else if( ldt.lt.nb ) then
     info = -10
  else if( ldc.lt.max( 1, m ) ) then
     info = -12
  end if

  if( info.ne.0 ) then
     call xerbla( 'dgemqrt', -info )
     return
  end if

  if( m.eq.0 .or. n.eq.0 .or. k.eq.0 ) return

  st = 1
  do while(stair(st)-ofs.le.0)
     st=st+1
     if(st.gt.k) return
  end do
  en = min(st+stair(k)-ofs-1,k)

  if( left .and. tran ) then

     i = 1
     do j = st, en, nb
        ib = min( nb, en-j+1 )
        im = max(min(stair(j+ib-1)-ofs,m) - i + 1,ib)
        if(im.gt.0) call _xlarfb( 'l', _qrm_transp, 'f', 'c', im, n, ib, &
             & v( i, j ), ldv, t( 1, j ), ldt, &
             & c( i, 1 ), ldc, work, ldwork )
        i = i+ib
     end do

  else if( right .and. notran ) then

     write(*,'("_qrm_gemqrt: not implemented")')
     info = -13

  else if( left .and. notran ) then

     kf = st+((en-st)/nb)*nb
     i  = kf-st+1
     do j = kf, st, -nb
        ib = min( nb, en-j+1 )
        im = max(min(stair(j+ib-1),m) - i + 1,ib)
        call _xlarfb( 'l', 'n', 'f', 'c', im, n, ib, &
             & v( i, j ), ldv, t( 1, j ), ldt, &
             & c( i, 1 ), ldc, work, ldwork )
        i = i-nb
     end do

  else if( right .and. tran ) then

     write(*,'("_qrm_gemqrt: not implemented")')
     info = -13

  end if

  return

end subroutine _qrm_gemqrt
