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
!> @file test_qrm.F90
!!
!! An example program which reads parameters from a file, a matrix from a MM file and
!! solves the problem. Run this program as
!!
!! ./_qrm_test < input.txt
!!
!! (where _ is d, s, c or z) using the provided input.txt file.
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version 0.0.1$
!! $Revision: 2244 $
!!
!!##############################################################################################

#include "qrm_common.h"

program _qrm_test

  use _qrm_mod
  implicit none

  type(_qrm_spmat_type)          :: qrm_mat

  integer                        :: INFO, nargs, i, nrhs, ounit, n
  character                      :: matfile*30='', transp 
  _qrm_data, allocatable, target :: b(:,:), x(:,:), r(:,:)
  integer, pointer               :: tmp(:)
  real(kind(1.d0))               :: t1, ta, tf, ts
  _qrm_real                      :: anrm
  _qrm_real, allocatable         :: rnrm(:), onrm(:), bnrm(:), xnrm(:)
  integer                        :: iseed(4)=(/1,1,1,1/)

  call qrm_init()
  
  ! initialize the matrix data structure. 
  call qrm_spmat_init(qrm_mat)

  ! read the input file
  call qrm_read_parms(qrm_mat, matfile, n, nrhs)
  call qrm_get('qrm_ounit', ounit)

  if(ounit.gt.0) write(ounit,'(30("="))')

  ! read the matrix. This subroutine has an overloaded interface
  ! and thus the type/precision of the data read will match that
  ! of qrm_mat
  if(matfile.eq.'Mitt3D_27') then
     call genmat_Mitt3D_27(qrm_mat, n)
  else
     call qrm_readmat(matfile, qrm_mat, .true., info)
     __QRM_INFO_CHECK(info,'qrm_test','qrm_readmat',10)
  end if

  write(ounit,'("Matrix ready. M:",i7,"  N:",i7,"  NZ:",i10)')qrm_mat%m, qrm_mat%n, qrm_mat%nz
  write(ounit,'(" ")')
  ! call _qrm_tikhonov(qrm_mat, 1d-8)
  
  if(qrm_mat%m .ge. qrm_mat%n) then
     transp='n'
  else
     if(ounit.gt.0) write(ounit,'("Transpose")')
     transp=_qrm_transp
  end if

  if(ounit.gt.0) write(ounit,'("Starting Analysis")')
  t1 = qrm_swtime()
  call qrm_analyse(qrm_mat, transp, info)
  ta = qrm_swtime()-t1
  __QRM_INFO_CHECK(info,'qrm_test','qrm_analyse',10)

  if(ounit.gt.0) write(ounit,'("  Estimated nonzeroes in R            : ",i20)')qrm_mat%gstats(qrm_e_nnz_r_)
  if(ounit.gt.0) write(ounit,'("  Estimated nonzeroes in H            : ",i20)')qrm_mat%gstats(qrm_e_nnz_h_)
  if(ounit.gt.0) write(ounit,'("  Estimated total flops at facto      : ",i20)')qrm_mat%gstats(qrm_e_facto_flops_)
  if(ounit.gt.0) write(ounit,'("  Estimated memory peak at facto (MB) : ",f20.3)')real(qrm_mat%gstats(qrm_e_facto_mempeak_))/1e6
  
  if(ounit.gt.0) write(ounit,'("Starting Factorization")')
  t1 = qrm_swtime()
  call qrm_factorize(qrm_mat, transp, info)
  tf = qrm_swtime()-t1
  __QRM_INFO_CHECK(info,'qrm_test','qrm_factorize',10)

  if(qrm_mat%icntl(qrm_keeph_) .eq. qrm_yes_) then
     if(ounit.gt.0) write(ounit,'("Starting Solve")')
     call qrm_alloc(b, qrm_mat%m, nrhs)
     call qrm_alloc(r, qrm_mat%m, nrhs)
     call qrm_alloc(x, qrm_mat%n, nrhs)
     
     call qrm_alloc(bnrm, nrhs)
     call qrm_alloc(rnrm, nrhs)
     call qrm_alloc(xnrm, nrhs)
     call qrm_alloc(onrm, nrhs)
     
     call _xlarnv(2, iseed, size(b), b(1,1))
     r = b

     t1 = qrm_swtime()
     if(transp .eq. 'n') then
        call qrm_apply(qrm_mat, _qrm_transp, b)
        call qrm_solve(qrm_mat, 'n', b, x)
     else if(transp .eq. _qrm_transp) then
        call qrm_solve(qrm_mat, _qrm_transp, b, x)
        call qrm_apply(qrm_mat, 'n', x)
     end if
     ts = qrm_swtime()-t1
     
     ! compute the residual
     call qrm_residual_norm(qrm_mat, r, x, rnrm)
     call qrm_vecnrm(x, size(x,1), '2', xnrm)
     call qrm_vecnrm(b, size(b,1), '2', bnrm)
     call qrm_matnrm(qrm_mat, 'f', anrm)
     call qrm_residual_orth(qrm_mat, r, onrm)   
     if(ounit.gt.0) then
        write(ounit,'(" ")')

        write(ounit,'("||A||          = ",e10.2)')anrm
        write(ounit,'(" ")')
        write(ounit,'("             ||b||         ||x||         ||r||/||A||   ||A^tr||/||r||")')
        write(ounit,'("---------------------------------------------------------------------")')
        do i=1, nrhs
           write(ounit,'("RHS ",i3,"  : ",4(e10.2,4x))')i,bnrm(i),xnrm(i),rnrm(i),onrm(i)
        end do
     end if
     
     call qrm_dealloc(b)
     call qrm_dealloc(r)
     call qrm_dealloc(x)
     call qrm_dealloc(bnrm)
     call qrm_dealloc(rnrm)
     call qrm_dealloc(xnrm)
     call qrm_dealloc(onrm)
  end if

10 continue

  call qrm_spmat_destroy(qrm_mat, all=.true.)

  if(ounit.gt.0) write(ounit,'(" ")')
  if(ounit.gt.0) write(ounit,'("Done.")')
  if(ounit.gt.0) write(ounit,'(" ")')
  if(ounit.gt.0) write(ounit,'(" ")')
 
  if(ounit.gt.0) write(ounit,'("  Time to do the analysis  : ",es10.3)')ta
  if(ounit.gt.0) write(ounit,'("  Time to do the facto     : ",es10.3)')tf
  if(ounit.gt.0) write(ounit,'("  Time to compute solution : ",es10.3)')ts
  if(ounit.gt.0) write(ounit,'("  Nonzeroes in R           : ",i20)')qrm_mat%gstats(qrm_nnz_r_)
  if(ounit.gt.0) write(ounit,'("  Nonzeroes in H           : ",i20)')qrm_mat%gstats(qrm_nnz_h_)
  if(ounit.gt.0) write(ounit,'("  Total unallocated memory : ",i20)')qrm_tot_mem
  if(ounit.gt.0) write(ounit,'("  Memory peak              : ",f9.3," MB")') &
       &real(qrm_max_mem,kind(1.d0))/1e6

  call qrm_finalize()
  stop

contains

  subroutine qrm_read_parms(qrm_mat, matfile, n, nrhs)
    ! Function: read_parms
    !
    ! *Input*:
    !
    ! *Output*:
    !
    use _qrm_spmat_mod
    implicit none

    type(_qrm_spmat_type)  :: qrm_mat
    character(len=*)       :: matfile
    integer                :: n, nrhs

    character              :: line*100, key*30, str*90
    integer                :: ival
    real(kind(1.d0))       :: rval

    nrhs = 1
    
    do
       read(*,'(a)')line
       read(line,*)key
       select case(key)
       case('end')
          exit
       case ('matfile')
          read(line,*)key, str
          str = adjustl(str)
          matfile = str(1:len_trim(str))
          if(matfile.eq.'Mitt3D_27') then
             read(line,*)key, str, n
          end if
       case ('qrm_ordering','qrm_keeph','qrm_nb','qrm_mb','qrm_bh',&
            & 'qrm_ib','qrm_rhsnb')
          read(line,*)key,ival
          call qrm_set(qrm_mat, key, ival)
       case ('qrm_mem_relax')
          read(line,*)key,rval
          call qrm_set(qrm_mat, key, rval)
       case ('qrm_ounit','qrm_eunit','qrm_dunit')
          read(line,*)key,ival
          call qrm_set(key, ival)
       case ('nrhs')
          read(line,*)key,nrhs
       case default
          if(ounit.gt.0) write(ounit,'("Unknown parameter ",a20)')key
       end select

    end do

    return
    
  end subroutine qrm_read_parms


  subroutine genmat_Mitt3D_27(qrm_mat, n)
    implicit none

    type(_qrm_spmat_type)  :: qrm_mat
    integer                :: n

    integer :: i, j, k, ig, idx, ounit
    integer, pointer :: irn(:), jcn(:)
    _qrm_data, pointer :: val(:)

    call qrm_get('qrm_ounit', ounit)

    qrm_mat%n  = (n+2)*(n+2)*(n+2);
    qrm_mat%m  = n*n*n;
    qrm_mat%nz = 27*n*n*n;

    if(ounit.gt.0)write(ounit,'("Generating Mittelmann Jacobian matrix for constraints of a")')
    if(ounit.gt.0)write(ounit,'("3D control problem with Dirichlet boundary conditions.")')
    ! if(ounit.gt.0)Write(ounit,'("N=",i3,"   ===>  m=",i7,"  n=",i7,"  nz=",i10)')n,qrm_mat%m, qrm_mat%n, qrm_mat%nz
    
    call qrm_alloc(qrm_mat%irn, qrm_mat%nz)
    call qrm_alloc(qrm_mat%jcn, qrm_mat%nz)
    call qrm_alloc(qrm_mat%val, qrm_mat%nz)
    irn => qrm_mat%irn
    jcn => qrm_mat%jcn
    val => qrm_mat%val

    idx = 1
    ig = 1
    do i=1, n
       do j=1, n
          do k=1, n

             irn(idx) = ig; jcn(idx) = y_index(n ,i-1 ,j   ,k  )  ; val(idx) = -16; idx=idx+1; ! y(i-1 ,j   ,k  )     
             irn(idx) = ig; jcn(idx) = y_index(n ,i   ,j   ,k  )  ; val(idx) = 200; idx=idx+1; ! y(i   ,j   ,k  )       
             irn(idx) = ig; jcn(idx) = y_index(n ,i+1 ,j   ,k  )  ; val(idx) = -16; idx=idx+1; ! y(i+1 ,j   ,k  )     
             irn(idx) = ig; jcn(idx) = y_index(n ,i-1 ,j-1 ,k  )  ; val(idx) = -8 ; idx=idx+1; ! y(i-1 ,j-1 ,k  )   
             irn(idx) = ig; jcn(idx) = y_index(n ,i   ,j-1 ,k  )  ; val(idx) = -16; idx=idx+1; ! y(i   ,j-1 ,k  )     
             irn(idx) = ig; jcn(idx) = y_index(n ,i+1 ,j-1 ,k  )  ; val(idx) = -8 ; idx=idx+1; ! y(i+1 ,j-1 ,k  )   
             irn(idx) = ig; jcn(idx) = y_index(n ,i-1 ,j+1 ,k  )  ; val(idx) = -8 ; idx=idx+1; ! y(i-1 ,j+1 ,k  )   
             irn(idx) = ig; jcn(idx) = y_index(n ,i   ,j+1 ,k  )  ; val(idx) = -16; idx=idx+1; ! y(i   ,j+1 ,k  )     
             irn(idx) = ig; jcn(idx) = y_index(n ,i+1 ,j+1 ,k  )  ; val(idx) = -8 ; idx=idx+1; ! y(i+1 ,j+1 ,k  )   
             irn(idx) = ig; jcn(idx) = y_index(n ,i-1 ,j   ,k-1)  ; val(idx) = -8 ; idx=idx+1; ! y(i-1 ,j   ,k-1)   
             irn(idx) = ig; jcn(idx) = y_index(n ,i   ,j   ,k-1)  ; val(idx) = -16; idx=idx+1; ! y(i   ,j   ,k-1)     
             irn(idx) = ig; jcn(idx) = y_index(n ,i+1 ,j   ,k-1)  ; val(idx) = -8 ; idx=idx+1; ! y(i+1 ,j   ,k-1)   
             irn(idx) = ig; jcn(idx) = y_index(n ,i-1 ,j-1 ,k-1)  ; val(idx) = -1 ; idx=idx+1; ! y(i-1 ,j-1 ,k-1) 
             irn(idx) = ig; jcn(idx) = y_index(n ,i   ,j-1 ,k-1)  ; val(idx) = -8 ; idx=idx+1; ! y(i   ,j-1 ,k-1)   
             irn(idx) = ig; jcn(idx) = y_index(n ,i+1 ,j-1 ,k-1)  ; val(idx) = -1 ; idx=idx+1; ! y(i+1 ,j-1 ,k-1) 
             irn(idx) = ig; jcn(idx) = y_index(n ,i-1 ,j+1 ,k-1)  ; val(idx) = -1 ; idx=idx+1; ! y(i-1 ,j+1 ,k-1) 
             irn(idx) = ig; jcn(idx) = y_index(n ,i   ,j+1 ,k-1)  ; val(idx) = -8 ; idx=idx+1; ! y(i   ,j+1 ,k-1)   
             irn(idx) = ig; jcn(idx) = y_index(n ,i+1 ,j+1 ,k-1)  ; val(idx) = -1 ; idx=idx+1; ! y(i+1 ,j+1 ,k-1) 
             irn(idx) = ig; jcn(idx) = y_index(n ,i-1 ,j   ,k+1)  ; val(idx) = -8 ; idx=idx+1; ! y(i-1 ,j   ,k+1)   
             irn(idx) = ig; jcn(idx) = y_index(n ,i   ,j   ,k+1)  ; val(idx) = -16; idx=idx+1; ! y(i   ,j   ,k+1)     
             irn(idx) = ig; jcn(idx) = y_index(n ,i+1 ,j   ,k+1)  ; val(idx) = -8 ; idx=idx+1; ! y(i+1 ,j   ,k+1)   
             irn(idx) = ig; jcn(idx) = y_index(n ,i-1 ,j-1 ,k+1)  ; val(idx) = -1 ; idx=idx+1; ! y(i-1 ,j-1 ,k+1) 
             irn(idx) = ig; jcn(idx) = y_index(n ,i   ,j-1 ,k+1)  ; val(idx) = -8 ; idx=idx+1; ! y(i   ,j-1 ,k+1)   
             irn(idx) = ig; jcn(idx) = y_index(n ,i+1 ,j-1 ,k+1)  ; val(idx) = -1 ; idx=idx+1; ! y(i+1 ,j-1 ,k+1) 
             irn(idx) = ig; jcn(idx) = y_index(n ,i-1 ,j+1 ,k+1)  ; val(idx) = -1 ; idx=idx+1; ! y(i-1 ,j+1 ,k+1) 
             irn(idx) = ig; jcn(idx) = y_index(n ,i   ,j+1 ,k+1)  ; val(idx) = -8 ; idx=idx+1; ! y(i   ,j+1 ,k+1)   
             irn(idx) = ig; jcn(idx) = y_index(n ,i+1 ,j+1 ,k+1)  ; val(idx) = -1 ; idx=idx+1; ! y(i+1 ,j+1 ,k+1) 

             ig=ig+1
             
          end do
       end do
    end do

    ! do i=1, qrm_mat%nz
       ! write(*,*)qrm_mat%irn(i),qrm_mat%jcn(i),qrm_mat%val(i)
    ! end do
    
    return
  end subroutine genmat_Mitt3D_27
  

  function y_index(n, i, j, k)
    integer :: n, i, j, k, y_index
    
    y_index = k + (n+2)*j + (n+2)*(n+2)*i + 1
    return
  end function y_index


  subroutine _qrm_tikhonov(qrm_mat, gamma)
    implicit none

    type(_qrm_spmat_type)  :: qrm_mat
    _qrm_real              :: gamma

    integer                :: i
    integer, parameter     :: ione=1
    _qrm_real              :: fnorma, _xnrm2

    fnorma = _xnrm2(qrm_mat%nz, qrm_mat%val(1), ione)
    call qrm_get('qrm_ounit', ounit)

    if(ounit.gt.0) write(ounit,'("Tikhonov regularization with gamma=",es10.2)')gamma*fnorma

    call qrm_realloc(qrm_mat%irn, qrm_mat%nz+min(qrm_mat%m,qrm_mat%n), copy=.true.)
    call qrm_realloc(qrm_mat%jcn, qrm_mat%nz+min(qrm_mat%m,qrm_mat%n), copy=.true.)
    call qrm_realloc(qrm_mat%val, qrm_mat%nz+min(qrm_mat%m,qrm_mat%n), copy=.true.)

    do i=1, min(qrm_mat%m,qrm_mat%n)
       qrm_mat%val(qrm_mat%nz+i) = gamma*fnorma
       if(qrm_mat%m .ge. qrm_mat%n) then
          qrm_mat%irn(qrm_mat%nz+i) = qrm_mat%m+i
          qrm_mat%jcn(qrm_mat%nz+i) = i
       else
          qrm_mat%jcn(qrm_mat%nz+i) = qrm_mat%n+i
          qrm_mat%irn(qrm_mat%nz+i) = i
       end if
    end do

    qrm_mat%nz = qrm_mat%nz+min(qrm_mat%m,qrm_mat%n)
    if(qrm_mat%m .ge. qrm_mat%n) then
       qrm_mat%m = qrm_mat%m + qrm_mat%n
    else
       qrm_mat%n = qrm_mat%n + qrm_mat%m
    end if
    return
  end subroutine _qrm_tikhonov

  
end program _qrm_test







