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

program _qrm_front

  use _qrm_mod
  use _qrm_methods_mod
  use qrm_string_mod
  implicit none

  type(_qrm_spmat_type), target  :: qrm_mat
  _qrm_data, allocatable         :: a(:,:), b(:), x(:), r(:), t(:,:), work(:)
  type(qrm_dscr_type)            :: qrm_dscr

  integer                        :: ierr, nargs, i, j, k, ounit, it, nth, pid
  type(_qrm_front_type), pointer :: front
  character                      :: transp
  integer                        :: m, n, mb, nb, ib, bh, warm
  real(kind(1.d0))               :: t1, ta, tf, ts, gflops
  character(LEN=20)              :: str
  integer                        :: iseed(4)
  _qrm_real                      :: rnrm, onrm, bnrm, xnrm, _xnrm2, err, anrm


  call get_command_argument(1,value=str,length=i)
  read(str(1:i),*)m; qrm_mat%m = m
  call get_command_argument(2,value=str,length=i)
  read(str(1:i),*)n; qrm_mat%n = n
  call get_command_argument(3,value=str,length=i)
  read(str(1:i),*)mb
  call get_command_argument(4,value=str,length=i)
  read(str(1:i),*)nb
  call get_command_argument(5,value=str,length=i)
  read(str(1:i),*)ib
  call get_command_argument(6,value=str,length=i)
  read(str(1:i),*)bh
  call get_command_argument(7,value=str,length=i)
  read(str(1:i),*)nth
  call get_command_argument(8,value=str,length=i)
  read(str(1:i),*)warm

  call qrm_init(nth)
#if defined(have_starpu) 
  call starpu_f_fxt_stop_profiling()
#endif
  
  ! initialize the control data structure.
  call qrm_spmat_init(qrm_mat)

  call qrm_set('qrm_ounit', -6)
  call qrm_set('qrm_eunit', -6)
  call qrm_set(qrm_mat, 'qrm_ordering', qrm_natural_)
  call qrm_set(qrm_mat, 'qrm_keeph', 1)
  call qrm_set(qrm_mat, 'qrm_mem_relax', -1.d0)
  call qrm_set(qrm_mat, 'qrm_nb', nb)
  call qrm_set(qrm_mat, 'qrm_ib', ib)
  call qrm_set(qrm_mat, 'qrm_mb', mb)
  call qrm_set(qrm_mat, 'qrm_bh', bh)


  call qrm_alloc(a, m, n)
  iseed = (/1,1,1,1/)
  call _xlarnv(2, iseed, m*n, a(1, 1))

  if(qrm_mat%m .ge. qrm_mat%n) then
     transp='n'
     gflops = (2.d0*(real(m)-real(n)/3.d0)*real(n)**2)/1e9
  else
     transp=_qrm_transp
     gflops = (2.d0*(real(n)-real(m)/3.d0)*real(m)**2)/1e9
  end if

  
  call fake_analysis(qrm_mat, transp)
  call fake_facto_init(qrm_mat, transp)

  front => qrm_mat%fdata%front_list(1)
  
  call _qrm_activate_front(qrm_mat, front)
  call _qrm_init_front(qrm_mat, front)

  call qrm_dscr_init(qrm_dscr)

  call fake_asm(front, a, transp)
  
#if defined(have_starpu) 
  call starpu_f_fxt_start_profiling()
#endif

  t1 = qrm_swtime()
  call _qrm_factorize_front(qrm_dscr, front, qrm_mat%fdata%work, 1)
  call qrm_barrier(qrm_dscr)
  tf = qrm_swtime()-t1
  
#if defined(have_starpu)
  call starpu_f_fxt_stop_profiling()
#endif

  call _qrm_clean_front(qrm_mat, front)
  call qrm_barrier()
  qrm_mat%fdata%ok = .true.
  
  call qrm_alloc(b, qrm_mat%m)
  call qrm_alloc(r, qrm_mat%m)
  call qrm_alloc(x, qrm_mat%n)
  
  call _xlarnv(2, iseed, m, b(1))
  r = b
  bnrm = _xnrm2(m, b, 1)

  if(transp .eq. 'n') then
     call qrm_apply(qrm_mat, _qrm_transp, b)
     call qrm_solve(qrm_mat, 'n', b, x)
  else if(transp .eq. _qrm_transp) then
     call qrm_solve(qrm_mat, _qrm_transp, b, x)
     call qrm_apply(qrm_mat, 'n', x)
  end if

  call qrm_barrier()

  call _xgemv('n', m, n, -_qrm_one, a, m, x, 1, _qrm_one, r, 1)
  rnrm = _xnrm2(m, r, 1)
  if(transp .eq. 'n') then
     call _xgemv(_qrm_transp, m, n, _qrm_one, a, m, r, 1, _qrm_zero, x, 1)
     onrm = _xnrm2(n, x, 1)
     err  = onrm/bnrm
  else
     err = rnrm
  end if
  
  write(*,'(i2, 3x, i5,3x,i5,3x,i3,3x,i3,3x,i3,3x,i3,3x)',advance='no')&
       & nth, m, n, &
       & qrm_mat%icntl(qrm_mb_), &
       & qrm_mat%icntl(qrm_nb_), &
       & qrm_mat%icntl(qrm_ib_), &
       & qrm_mat%icntl(qrm_bh_)
  write(*,'(f10.2,3x)',advance='no')gflops
  write(*,'(es10.3,2x)',advance='no')tf
  write(*,'(f6.2,2x)',advance='no')gflops/tf
  write(*,'(es10.3,2x)')err

  call qrm_spmat_destroy(qrm_mat, all=.true.)

  stop

contains

  subroutine fake_analysis(qrm_mat, transp)
    use _qrm_fdata_mod
    implicit none

    type(_qrm_spmat_type), target   :: qrm_mat
    character                       :: transp
    type(qrm_adata_type), pointer   :: adata

    integer :: i, j, m, n

    call qrm_adata_init(qrm_mat%adata)
    adata => qrm_mat%adata

    if(transp.eq._qrm_transp) then
       m = qrm_mat%n
       n = qrm_mat%m
    else
       m = qrm_mat%m
       n = qrm_mat%n
    end if
    
    call qrm_alloc(adata%cperm,    n)
    call qrm_alloc(adata%rperm,    m)
    call qrm_alloc(adata%cp_ptr,   3)
    call qrm_alloc(adata%rc,       2)
    call qrm_alloc(adata%parent,   2)
    call qrm_alloc(adata%childptr, 3)
    call qrm_alloc(adata%child,    1)
    call qrm_alloc(adata%nfrows,   2)
    call qrm_alloc(adata%stair,    2)
    call qrm_alloc(adata%small,    2)
    call qrm_alloc(adata%fcol,     n)
    call qrm_alloc(adata%fcol_ptr, 3)
    call qrm_alloc(adata%leaves,   1)
    call qrm_alloc(adata%torder,   2)
    call qrm_alloc(adata%asize,    2)
    call qrm_alloc(adata%csize,    2)

    adata%cperm       = (/(j, j=1, n)/)
    adata%rperm       = (/(i, i=1, m)/)

    adata%cp_ptr(1)   = 1; adata%cp_ptr(2) = n+1; adata%cp_ptr(3)=n+1

    adata%rc(1)       = n; adata%rc(2)     = 0
    adata%nfrows(1)   = m; adata%nfrows(2) = 0

    adata%parent(1)   = 2; adata%parent(2) = 0

    adata%childptr(1) = 1; adata%childptr(2) = 1; adata%childptr(3) = 2
    adata%child(1)    = 1
    
    adata%stair(1)    = m; adata%stair(2)    = m

    adata%small       = 0

    adata%fcol        = (/(j, j=1, n)/)
    adata%fcol_ptr(1) = 1; adata%fcol_ptr(2) = n+1; adata%fcol_ptr(3) = n+1

    adata%leaves(1)   = 1

    adata%torder      = (/1, 2/)

    adata%asize       = 1
    adata%csize       = 1

    adata%nnodes      = 2
    adata%ok          = .true.

    return
  end subroutine fake_analysis

  subroutine fake_facto_init(qrm_mat, transp)
    use _qrm_fdata_mod
    implicit none
    
    type(_qrm_spmat_type), target   :: qrm_mat
    character                       :: transp
    type(_qrm_fdata_type), pointer   :: fdata
    type(qrm_adata_type), pointer   :: adata
    type(_qrm_front_type), pointer :: front
    
    integer :: i, j, node
    
    call qrm_fdata_init(qrm_mat%fdata)
    fdata => qrm_mat%fdata
    adata => qrm_mat%adata

    fdata%nfronts = adata%nnodes

    
    allocate(fdata%front_list(adata%nnodes))

    call qrm_facto_mem_init(fdata%ma, huge(1_8))

    
#if defined(have_starpu)
    call starpu_f_matrix_data_register(fdata%work%hdl, -1, c_null_ptr, &
         & qrm_mat%icntl(qrm_nb_), &
         & qrm_mat%icntl(qrm_nb_), &
         & qrm_mat%icntl(qrm_nb_), &
         & int(_qrm_sizeof_data,kind=c_size_t))
#else
    allocate(fdata%work%c(qrm_mat%icntl(qrm_nb_),qrm_mat%icntl(qrm_nb_)))
#endif

    front => fdata%front_list(1)
#if defined(have_starpu)
    call starpu_f_void_data_register(front%sym_handle)
#endif
    front%anrows = adata%stair(1)
    call qrm_alloc(front%aiptr, front%anrows+1)
    call qrm_alloc(front%ajcn,  front%anrows+adata%rc(1)-1)
    call qrm_alloc(front%aval,  front%anrows+adata%rc(1)-1)
    front%aiptr(1:front%anrows)                         = (/(j, j=1, front%anrows)/)
    front%aiptr(front%anrows+1)                         = front%anrows+adata%rc(1)
    front%ajcn(1:front%anrows-1)                        = 1
    front%ajcn(front%anrows:front%anrows+adata%rc(1)-1) = (/(j, j=1, adata%rc(1))/)
    front%aval                                          = _qrm_zero
    front%num = 1
    
    front => fdata%front_list(2)
#if defined(have_starpu)
    call starpu_f_void_data_register(front%sym_handle)
#endif
    front%anrows = 0
    front%num    = 2
    front%m      = 0
    front%n      = 0
    front%ne     = 0
    
    qrm_mat%gstats(qrm_nnz_r_) = 0
    qrm_mat%gstats(qrm_nnz_h_) = 0
    
    return
  end subroutine fake_facto_init

  
  subroutine fake_asm(front, a, transp)
    implicit  none
    
    type(_qrm_front_type) :: front
    _qrm_data             :: a(:,:)
    character             :: transp

    integer :: i, j, mm, nn

    if((transp.eq._qrm_transp) .or. (transp.eq.'c')) then
       colst: do j=1, front%nc
          nn = min(front%nb,front%n-(j-1)*front%nb) 
          rowst: do i=1, front%nr
             mm = min(front%stair(min(j*front%nb,front%n)) - (i-1)*front%mb, front%mb)
             if(mm .le. 0) cycle colst
             call _xcopy_transpose(transp, front%bc(i,j)%c, &
                  & a((j-1)*front%nb+1:(j-1)*front%nb+nn, (i-1)*front%mb+1:(i-1)*front%mb+mm))
          end do rowst
       end do colst
    else
       cols: do j=1, front%nc
          nn = min(front%nb,front%n-(j-1)*front%nb) 
          rows: do i=1, front%nr
             mm = min(front%stair(min(j*front%nb,front%n)) - (i-1)*front%mb, front%mb)
             if(mm .le. 0) cycle cols
             ! front%bc(i,j)%c = a((i-1)*front%mb+1:(i-1)*front%mb+mm,(j-1)*front%nb+1:(j-1)*front%nb+nn)
             call _xcopy_transpose(transp, front%bc(i,j)%c, &
                  & a((i-1)*front%mb+1:(i-1)*front%mb+mm,(j-1)*front%nb+1:(j-1)*front%nb+nn))
          end do rows
       end do cols
    end if
    return
  end subroutine fake_asm

  subroutine _xcopy_transpose(transp, a, b)
    implicit none
    character             :: transp
    _qrm_data             :: a(:,:), b(:,:)

    integer :: i, j

    if(transp.eq._qrm_transp) then
       do i=1, size(a,1)
          do j=1, size(a,2)
#if defined(zprec) || defined(cprec)
             a(i,j) = conjg(b(j,i))
#else
             a(i,j) = b(j,i)
#endif
          end do
       end do
    else
       do i=1, size(a,1)
          do j=1, size(a,2)
             a(i,j) = b(i,j)
          end do
       end do
    end if
    
    return
  end subroutine _xcopy_transpose

end program _qrm_front

