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
!> @file qrm_utils_mod.F90
!! This file contains a module with generic interfaces for a number of auxiliary tools
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


!> @brief This module contains generic interfaces for a number of auxiliary tools
module _qrm_utils_mod
  use iso_c_binding

  !> @brief Generic interface for the @link ::_qrm_readmat @endlink
  !> routine
  interface qrm_readmat
     subroutine _qrm_readmat(matfile, qrm_mat, fakec, info)
       use _qrm_spmat_mod
       character, intent(in)                :: matfile*30
       type(_qrm_spmat_type), intent(inout) :: qrm_mat
       logical, optional                    :: fakec
       integer, optional                    :: info
     end subroutine _qrm_readmat
  end interface qrm_readmat

  !> @brief Generic interface for the @link ::_qrm_matmul2d @endlink
  !> and @link ::_qrm_matmul1d @endlink routines
  interface qrm_matmul
     procedure :: _qrm_matmul2d, _qrm_matmul1d
  end interface qrm_matmul

  !> @brief Generic interface for the @link ::_qrm_matmul2d @endlink
  !> and @link ::_qrm_matmul1d @endlink routines
  interface _qrm_matmul
     subroutine _qrm_matmul2d(qrm_mat, transp, alpha, x, beta, y)
       use _qrm_spmat_mod
       type(_qrm_spmat_type), target :: qrm_mat
       _qrm_data, intent(out)        :: y(:,:)
       _qrm_data, intent(in)         :: x(:,:)
       _qrm_data, intent(in)         :: alpha, beta
       character(len=*)              :: transp
     end subroutine _qrm_matmul2d
     subroutine _qrm_matmul1d(qrm_mat, transp, alpha, x, beta, y)
       use _qrm_spmat_mod
       type(_qrm_spmat_type)         :: qrm_mat
       _qrm_data, intent(out)        :: y(:)
       _qrm_data, intent(in)         :: x(:)
       _qrm_data, intent(in)         :: alpha, beta
       character(len=*)              :: transp
     end subroutine _qrm_matmul1d
  end interface _qrm_matmul


  !> @brief Generic interface for the @link ::_qrm_vecnrm2d @endlink
  !> and @link ::_qrm_vecnrm1d @endlink routines
  interface qrm_vecnrm
     procedure :: _qrm_vecnrm2d, _qrm_vecnrm1d
  end interface qrm_vecnrm

  !> @brief Generic interface for the @link ::_qrm_vecnrm2d @endlink
  !> and @link ::_qrm_vecnrm1d @endlink routines
  interface _qrm_vecnrm
     subroutine _qrm_vecnrm2d(vec, n, ntype, nrm, info)
       _qrm_data, intent(in)  :: vec(:,:)
       _qrm_real              :: nrm(:)
       integer, intent(in)    :: n
       character              :: ntype
       integer, optional      :: info
     end subroutine _qrm_vecnrm2d
     subroutine _qrm_vecnrm1d(vec, n, ntype, nrm, info)
       _qrm_data, intent(in)  :: vec(:)
       _qrm_real              :: nrm
       integer, intent(in)    :: n
       character              :: ntype
       integer, optional      :: info
     end subroutine _qrm_vecnrm1d
  end interface _qrm_vecnrm

  !> @brief Generic interface for the @link ::_qrm_remap_pnt @endlink
  !> routine
  interface qrm_remap_pnt
     subroutine _qrm_remap_pnt(arr1d, pnt2d, n)
       integer :: n
       _qrm_data, target  :: arr1d(1:n)
       _qrm_data, pointer :: pnt2d(:,:)
     end subroutine _qrm_remap_pnt
  end interface


  
  interface qrm_matnrm
     subroutine _qrm_matnrm(qrm_mat, ntype, nrm, info)
       use _qrm_spmat_mod
       type(_qrm_spmat_type), intent(in) :: qrm_mat
       _qrm_real                         :: nrm
       character                         :: ntype
       integer, optional                 :: info
     end subroutine _qrm_matnrm
  end interface qrm_matnrm


  interface qrm_get_front_mem
     module procedure _qrm_get_front_mem
  end interface qrm_get_front_mem
  
  interface qrm_get_front_flops
     module procedure _qrm_get_front_flops
  end interface qrm_get_front_flops
  

contains

  subroutine _qrm_get_front_mem(m, n, mb, nb, ib, bh, np, amem, cmem, hsize, rsize, storer, storeh, stair)
    use qrm_adata_mod
    use qrm_common_mod
    implicit none
    integer              :: m, n, mb, nb, ib, np, bh
    ! if stair is not present we assume that the node is small
    integer, optional    :: stair(:)
    integer(kind=8)      :: amem, cmem, hsize, rsize
    logical              :: storer, storeh
    
    integer(kind=8)      :: ne, cm, cn, cne, cnp, c, ibh, iib, imb, inb
    integer(kind=8)      :: nn, mm, i, j, k, nc, nr

    if(mb.lt.0) then
       imb = m
    else
       imb = min(mb,m)
    end if
    ne = min(m,n)
    inb = min(nb,n)
    iib = min(int(ib,kind=8),int(ne,kind=8))
    ibh = bh
        
    hsize = 0
    rsize = 0
    amem  = 0
    cmem  = 0

    if(min(m,n).le.0) return

    nc = (n-1)/inb + 1
    nr = (m-1)/imb + 1
    if(ibh.le.0) ibh = nr
    
    ! count only real/complex data for the moment
    ! copy/paste the allocation loop in qrm_activate_front

    cols: do j=1, nc
       nn = min(inb,int(n-(j-1)*inb,kind=8))

       if(present(stair)) then
          if((j-1)*inb .lt. ne) hsize = hsize + nn*(stair((j-1)*inb+nn)-(j-1)*inb)
       else
          if((j-1)*inb .lt. ne) hsize = hsize + nn*(m-(j-1)*inb)
       end if
       rsize = rsize + nn*min(int(np,kind=8),(j-1)*inb+nn)

       rows: do i=1, nr
          if(present(stair)) then
             mm = min(int(stair(min(j*inb,int(n,kind=8))) - (i-1)*imb,kind=8), imb)
          else
             mm = min(m - (i-1)*imb, imb)
          end if
          if(mm .le. 0) cycle cols
          amem = amem + mm*nn
          if((i-1)*imb+mm .ge. (j-1)*inb+1) then
             k = ((j-1)*inb)/imb+1
             if(mod(i-k,ibh).eq.0) then
                ! For those blocks where _geqrt is done, after the
                ! _geqrt operation the v vectors are copied into the
                ! front%v block along with the T matrices. This trick
                ! allows us to remove a false dependency between _gemqrt
                ! and _tpqrt tasks. Padding is added to avoid TRSMs on
                ! the GPU and thus the T matrices are literally stored
                ! on top of the V blocks:
                !   |TTTT|
                !   |V000|
                !   |VV00|
                !   |VVV0|
                !   |VVVV|
                hsize = hsize + iib*nn
                if(.not.present(stair)) then
                   mm = iib
                else
                   mm = mm - max((j-1)*inb +1 - (i-1)*imb,1_8)+1 + iib
                end if
                ! write(*,*)m,n,imb,inb,iib,mm,nn
                amem = amem + mm*nn
             end if
             if(present(stair))  then
                hsize = hsize + iib*nn
                amem = amem + iib*nn
             end if
          end if
       end do rows
    end do cols

    if(storeh.and.storer) then
       cmem = amem - (hsize+rsize)
    else if (storer) then
       cmem = amem - rsize
    else
       cmem = amem
    end if
    
    ! add one temp block for the clean
    amem = amem+2*inb*imb
    ! add one temp block for the clean
    cmem = cmem+2*inb*imb

    ! hsize and rsize are expressed in # of entries whereas
    ! amem and cmem are expressed in Bytes
    amem = amem*_qrm_sizeof_data
    cmem = cmem*_qrm_sizeof_data

    ! cols
    amem = amem+n*_qrm_sizeof_i
    ! rows
    amem = amem+m*_qrm_sizeof_i
    ! stair
    amem = amem+(n+1)*_qrm_sizeof_i
    ! colmap
    amem = amem+(n-np)*_qrm_sizeof_i
    ! rowmap
    amem = amem+(ne-np)*_qrm_sizeof_i

    if(.not.storer) then
       ! cols
       cmem = cmem+n*_qrm_sizeof_i
       ! rows
       cmem = cmem+m*_qrm_sizeof_i
       ! stair
       cmem = cmem+(n+1)*_qrm_sizeof_i
       ! colmap
       cmem = cmem+(n-np)*_qrm_sizeof_i
       ! rowmap
       cmem = cmem+(ne-np)*_qrm_sizeof_i
    end if
    
    ! the amount of memory freed by the clean routine may be
    ! negaitve. This is due to the fact that there may be some overlap
    ! between the H and R blocks. This is just a quick fix.
    if(cmem.lt.0) then
       amem = amem-cmem
       cmem = 0
    end if

    return
  end subroutine _qrm_get_front_mem



  subroutine _qrm_get_front_flops(m, n, stair, mb, nb, ib, flops)
    use qrm_common_mod
    implicit none
    integer :: m, n, mb, nb, ib
    integer :: stair(:)
    integer(kind=8) :: flops

    integer :: j, ne
    integer :: nn, mm, kk

    ne = min(m,n)
    flops = 0

    if(min(m,n).le.0) return
        
    do j=1, ne, ib
       kk = min(ib,ne-j+1)
       mm = stair(j+kk-1) - j+1
       flops = flops + qrm_count_realflops(mm, kk, kk, 'panel')
       nn = n - (j+kk) + 1
       if(nn.gt.0) flops = flops + qrm_count_realflops(mm, nn, kk, 'update')
    end do

  end subroutine _qrm_get_front_flops


  function qrm_compute_task_flops(optype, m, n, k, l, stair, ofsa, ofsb)
    implicit none

    integer         :: m, n, k, l, ofsa, ofsb, stair(:)
    character       :: optype*(*)
    integer(kind=8) :: qrm_compute_task_flops
    
    integer         :: i, j, im, in
    integer(kind=8) :: flops

    flops = 0

    select case(optype)
    case ('geqrt')
       i = 0
       do j=1, min(n,min(stair(n)-ofsa,m))
          im = min(max(stair(j)-ofsa,0),m)-i
          if(im.le.0) cycle
          in = n-j
          flops = flops + im*(3+4*in)
          i = i+1
       end do
    case ('gemqrt')
       i = 0
       do j=1, min(k,min(stair(k)-ofsa,m))
          im = min(max(stair(j)-ofsa,0),m)-i
          if(im.le.0) cycle
          flops = flops + im*4*n
          i = i+1
       end do
    case ('tpqrt')
       im = 0
       do j=1, n
          if(stair(j).lt.ofsb) cycle
          if(l.eq.m) then
             im = min(im+1,m)
          else if(l.eq.0) then
             im = min(m, stair(j)-ofsb)
          end if
          in = n-j
          flops = flops + (im+1)*(3+4*in)
       end do
    case ('tpmqrt')
       im = 0
       do j=1, k
          if(stair(j).lt.ofsb) cycle
          if(l.eq.m) then
             im = min(im+1,m)
          else if(l.eq.0) then
             im = min(m, stair(j)-ofsb)
          end if
          flops = flops + (im+1)*4*n
       end do
    end select

    if(flops.lt.0) write(*,*)'Error in flopcount'

    qrm_compute_task_flops = flops
    
    return
  end function qrm_compute_task_flops


end module _qrm_utils_mod







