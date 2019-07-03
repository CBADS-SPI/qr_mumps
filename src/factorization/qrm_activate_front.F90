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


#include "qrm_common.h"

subroutine _qrm_activate_front(qrm_mat, front, work, info)
  use _qrm_spmat_mod
  use _qrm_fdata_mod
  use qrm_mem_mod
  use _qrm_utils_mod
  use qrm_common_mod
#if defined (have_starpu)
  use starpu_f_mod
#endif
  implicit none

  type(_qrm_spmat_type), target  :: qrm_mat
  type(_qrm_front_type)          :: front
  integer, optional              :: work(:)
  integer, optional              :: info

  type(qrm_adata_type), pointer  :: adata
  type(_qrm_fdata_type), pointer :: fdata
  type(_qrm_front_type), pointer :: cfront
  integer                        :: i, j, k
  integer                        :: m, n, mb, nb, ib, ne, npiv
  integer                        :: cbc, cbr, fnum, father, nc, nr, mm, nn, bh, tmpmem

  ! error management
  integer                        :: err
  character(len=*), parameter    :: name='qrm_activate_front'


  err = 0

  ! tmpmem = qrm_tot_mem
  
  adata      => qrm_mat%adata
  fdata      => qrm_mat%fdata
  fnum       =  front%num
  front%m    =  adata%nfrows(fnum)
  front%n    =  adata%rc(fnum)
  ! The number of eliminations to be performed on the front f
  front%ne   = min(front%m,front%n)

  if( (front%n .le. 0) .or. (front%m .le. 0)) then
     ! nothing to do here. Mark the front as done
     front%np = 0
     front%nc = 0
     goto 9999
  end if

  front%ib   =  min(qrm_mat%icntl(qrm_ib_),front%ne)
  if(adata%small(fnum).eq.0) then
     front%nb   =  min(qrm_mat%icntl(qrm_nb_),front%n)
     if(qrm_mat%icntl(qrm_mb_).le.0) then
        front%mb   = front%m
     else
        front%mb   =  min(qrm_mat%icntl(qrm_mb_),front%m)
     end if
  else
     ! no blocking is done on small fronts because they are treated
     ! sequentially
     front%nb   =  min(qrm_mat%icntl(qrm_nb_),front%n)
     ! front%nb   =  front%n
     front%mb   =  front%m
  end if


  father     =  adata%parent(fnum)

  m          = front%m
  n          = front%n
  nb         = front%nb
  mb         = front%mb
  ib         = front%ib
  ne         = front%ne
  nc         = (n-1)/nb + 1
  nr         = (m-1)/mb + 1
  npiv       = min(adata%cp_ptr(fnum+1)-adata%cp_ptr(fnum), ne)

  ! The number of panel columns in f
  front%np   = max((front%ne-1)/nb + 1,0)
  ! The number of block-rows in f
  front%nr   = nr
  ! The number of block-columns in f
  front%nc   = nc
  ! The number of pivots in f
  front%npiv = min(adata%cp_ptr(fnum+1)-adata%cp_ptr(fnum),front%ne)

  if(qrm_mat%icntl(qrm_bh_).le.0) then
     front%bh = front%nr
  else
     front%bh = qrm_mat%icntl(qrm_bh_)
  end if

  bh         = front%bh
  
  cbr        = ne-npiv   ! The number of rows in CB
  cbc        = n-npiv    ! The number of cols in the cb

  if(err.eq.0) call qrm_alloc(front%cols,  n,   err)
  if(err.eq.0) call qrm_alloc(front%stair, n+1, err)
  if(err.eq.0) call qrm_alloc(front%rows,  m,   err)
  __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)

  ! initialize the column indices and build the colmap
  front%cols(1:front%n) = adata%fcol(adata%fcol_ptr(fnum): &
       & adata%fcol_ptr(fnum+1)-1)

  call _qrm_build_colmapping(adata, fdata, front, work, err)
  __QRM_INFO_CHECK(err, name,'qrm_build_colmapping',9999)

  call qrm_alloc(front%arowmap, front%anrows, err)
  __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)

  call _qrm_build_stair_rowmapping(adata, fdata, front, err)
  __QRM_INFO_CHECK(err, name,'qrm_build_stair_rowmapping',9999)

  ! allocate block-columns
  allocate(front%bc(nr,nc))
  allocate(front%t (nr,nc))
  allocate(front%t2(nr,nc))
  cols: do j=1, nc
     nn = min(nb,n-(j-1)*nb)
     rows: do i=1, nr
        mm = min(front%stair(min(j*nb,n)) - (i-1)*mb, mb)
        if(mm .le. 0) cycle cols
        call qrm_alloc( front%bc(i,j)%c, mm, nn, err )
        if(err.ne.0) exit cols
        if((i-1)*front%mb+mm .ge. (j-1)*front%nb+1) then
           k = ((j-1)*nb)/mb+1
           if(mod(i-k,bh).eq.0) then
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
              if(adata%small(fnum).eq.0) then
                 mm = mm - max((j-1)*nb +1 - (i-1)*mb,1)+1 + ib
              else
                 mm = ib
              end if
              call qrm_alloc( front%t(i,j)%c, mm, nn, err )
              if(err.ne.0) exit cols
           end if
           if(adata%small(fnum).eq.0) call qrm_alloc( front%t2(i,j)%c, ib, nn, err )
           if(err.ne.0) exit cols
        end if
     end do rows
  end do cols
  __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)

#if defined(have_starpu)
  if(adata%small(fnum).eq.0) call _qrm_front_data_register(front)
#endif

  ! write(*,'(i5,2x,i5," activ -- ",i10,2x,i10,2x,f7.2,2x,i10,2x,i10,2x,i10)')front%num, qrm_mat%adata%small(front%num), &
       ! & qrm_tot_mem-tmpmem, &
       ! & adata%asize(fnum), &
       ! & 100.d0*real(adata%asize(fnum)-(qrm_tot_mem-tmpmem))/real(qrm_tot_mem-tmpmem), &
       ! & qrm_tot_mem, qrm_max_mem, qrm_mat%fdata%ma%peak

9999 continue 
  ! should do cleanup
  if(present(info)) info = err
  return


contains

#if defined(have_starpu)
  subroutine _qrm_front_data_register(front)
    use _qrm_fdata_mod
    implicit none

    type(_qrm_front_type) :: front
 
    integer               :: mm, i, j, k

    front%bc(:,:)%hdl  = c_null_ptr
    front%t(:,:)%hdl   = c_null_ptr
    front%t2(:,:)%hdl  = c_null_ptr
    
    cols: do j=1, front%nc
       rows: do i=1, front%nr
          mm = min(front%stair(min(j*front%nb,front%n)) - (i-1)*front%mb, front%mb)
          if(mm .le. 0) cycle cols
          
          call starpu_f_matrix_data_register(front%bc(i,j)%hdl, 0, &
               & front%bc(i,j)%c, &
               & size(front%bc(i,j)%c,1), &
               & size(front%bc(i,j)%c,1), &
               & size(front%bc(i,j)%c,2))

          ! this is not needed for all the blocks in the front
          if((i-1)*front%mb+mm .ge. (j-1)*front%nb+1) then

             k = ((j-1)*nb)/mb+1
             if(mod(i-k,bh).eq.0) then
                call starpu_f_matrix_data_register(front%t(i,j)%hdl, 0, &
                     & front%t(i,j)%c, &
                     & size(front%t(i,j)%c,1), &
                     & size(front%t(i,j)%c,1), &
                     & size(front%t(i,j)%c,2))
             end if
             
             call starpu_f_matrix_data_register(front%t2(i,j)%hdl, 0, &
                  & front%t2(i,j)%c, &
                  & size(front%t2(i,j)%c,1), &
                  & size(front%t2(i,j)%c,1), &
                  & size(front%t2(i,j)%c,2))
          end if
       end do rows
    end do cols

  end subroutine _qrm_front_data_register
#endif

  subroutine _qrm_build_colmapping(adata, fdata, front, work, info)

    implicit none
    type(qrm_adata_type), target   :: adata
    type(_qrm_fdata_type), target  :: fdata
    type(_qrm_front_type)          :: front
    integer, target, optional      :: work(:)
    integer, optional              :: info

    integer                        :: nch, f, c, i, j, k, cn, cnpiv, ce
    integer, pointer               :: ccols(:), gcolmap(:)
    integer, allocatable           :: ptrs(:)
    type(_qrm_front_type), pointer :: cfront
    logical                        :: goodmem
    ! error management
    integer                        :: err
    character(len=*), parameter    :: name='qrm_activate_front'

    err = 0
    
    if(present(work)) then
       if(size(work).lt.qrm_mat%n) then
          goodmem = .true.
       else
          goodmem = .false.
          gcolmap => work
       end if
    else
       goodmem = .true.
    end if

    f = front%num
    nch = adata%childptr(f+1)-adata%childptr(f)

    if (.not.goodmem) then
    
       do j=1, front%n
          k = front%cols(j)
          gcolmap(k)=j
       end do
    
       do i=1, nch
          c = adata%child(adata%childptr(f)+i-1)
          cfront => fdata%front_list(c)
          cfront%n    = adata%rc(c)
          cfront%ne   = min(cfront%n,adata%nfrows(c))
          cfront%npiv = min(adata%cp_ptr(c+1)-adata%cp_ptr(c), cfront%ne)
         ! allocate(cfront%colmap(cfront%n-cfront%npiv))
          call qrm_alloc(cfront%colmap, cfront%n-cfront%npiv, err)
          __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)
          ccols => adata%fcol(adata%fcol_ptr(c):adata%fcol_ptr(c+1)-1)
          do j=cfront%npiv+1, cfront%n
             ! this is the column mapping on the child. cfront%colmap(k)=j
             ! means that the k-th column of cfront will be assembled into the
             ! j-th column of front
             cfront%colmap(j-cfront%npiv) = gcolmap(ccols(j))
          end do
       end do
       
       do i=1, front%anrows
          do j= front%aiptr(i), front%aiptr(i+1)-1
             front%ajcn(j) = gcolmap(front%ajcn(j))
          end do
       end do
       
    else

       do i=1, nch
          c = adata%child(adata%childptr(f)+i-1)
          cfront => fdata%front_list(c)
          cfront%n    = adata%rc(c)
          cfront%ne   = min(cfront%n,adata%nfrows(c))
          cfront%npiv = min(adata%cp_ptr(c+1)-adata%cp_ptr(c), cfront%ne)
          ! allocate(cfront%colmap(cfront%n-cfront%npiv))
          call qrm_alloc(cfront%colmap, cfront%n-cfront%npiv, err)
          __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)
       end do
       
       allocate(ptrs(max(nch,front%anrows)))
       
       ptrs = 1
       ! compute the mapping for the child nodes.  cfront%colmap(i)=j
       ! means that column i+cfront%npiv of the child has to be assembled
       ! into column j of the father.
       ! The instructions below rely on the fact that both front%cols and
       ! cfront%cols are sorted with respect to the global column
       ! permutation
       do j=1, front%n
          do i=1, nch
             c = adata%child(adata%childptr(f)+i-1)
             cfront => fdata%front_list(c)
             if(ptrs(i).gt.cfront%n-cfront%npiv) cycle
             
             ccols => adata%fcol(adata%fcol_ptr(c):adata%fcol_ptr(c+1)-1)
             
             if(ccols(cfront%npiv+ptrs(i)).eq.front%cols(j)) then
                cfront%colmap(ptrs(i)) = j
                ptrs(i) = ptrs(i)+1
             end if
             
          end do
       end do
       
       ! As for the coefficients of the original sparse matrix, we simply
       ! replace the corresponding global column indices with the column
       ! indices of the front. These instructions also rely on the
       ! assumption that front%ajcn is sorted wrt the global column
       ! permutation (which is doen in the qrm_factorization_init routine
       ptrs(1:front%anrows) = front%aiptr(1:front%anrows)
       do j=1, front%n
          do i=1, front%anrows
             if(ptrs(i).ge.front%aiptr(i+1)) cycle
             ! this while loop is necessary to cope with the presence of
             ! duplicates
             do while(front%ajcn(ptrs(i)).eq.front%cols(j))
                front%ajcn(ptrs(i)) = j
                ptrs(i) = ptrs(i)+1
                if(ptrs(i).ge.front%aiptr(i+1)) exit
             end do
          end do
       end do

       deallocate(ptrs)
       
    end if

9999 continue
    if(present(info)) info = err
    return

  end subroutine _qrm_build_colmapping




  subroutine _qrm_build_stair_rowmapping(adata, fdata, front, info)
    implicit none

    type(_qrm_front_type)          :: front
    type(qrm_adata_type)           :: adata
    type(_qrm_fdata_type), target  :: fdata
    integer, optional              :: info
    
    integer, allocatable           :: first(:)
    integer                        :: i, ne, npiv, f, c, p, row, roff

    ! error management
    integer                        :: err
    character(len=*), parameter    :: name='qrm_build_stair_rowmapping'

    err = 0
    
    if (front%num .eq. 1) then
       roff = 0
    else
       roff = adata%stair(front%num-1)
    end if
    
    if(front%ne .lt. front%ib) then
       front%stair = front%m
       ! build rowmap for the rows from the original matrix
       row = 0
       do i=1, front%anrows
          row = row+1
          front%arowmap(i) = row
          front%rows(row) = adata%rperm(roff+i)
       end do
       
       ! build rowmap for the rows from the children
       do p = adata%childptr(fnum), adata%childptr(fnum+1)-1
          c = adata%child(p)
          cfront => fdata%front_list(c)
          ! allocate(cfront%rowmap(cfront%ne-cfront%npiv))
          call qrm_alloc(cfront%rowmap, cfront%ne-cfront%npiv, err)
          __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)
               
          ! ne is the number of Householder vectors computed on the
          ! child c. npiv is the number of fully assembled pivots in c
          ne   = cfront%ne
          npiv = cfront%npiv
          
          ! this is the row mapping on the child. cfront%rowmap(k)=i
          ! means that the k-th row of cfront will be assembled into the
          ! i-th row of front
          do i=npiv+1, ne
             row = row+1
             cfront%rowmap(i-npiv) = row
          end do
       end do
    else
       front%stair = 0
       allocate(first(front%anrows))
       
       first = front%n+1
       ! count in the rows from the original matrix
       do i=1, front%anrows
          ! This instruction below relies on the assumption that in each
          ! row, the coefficients are sorted in increasing order wrt the
          ! global permutation
          first(i) = front%ajcn(front%aiptr(i))
          front%stair(first(i)+1) = front%stair(first(i)+1)+1
       end do

       ! count in rows coming from CBs
       do p = adata%childptr(fnum), adata%childptr(fnum+1)-1
          c = adata%child(p)
          cfront => fdata%front_list(c)
          
          ! ne is the number of Householder vectors computed on the
          ! child c. npiv is the number of fully assembled pivots in c
          ne = cfront%ne
          npiv = cfront%npiv
          
          ! count in all the rows on the CB of c
          do i=1, ne-npiv
             f = cfront%colmap(i)
             front%stair(f+1) = front%stair(f+1)+1
          end do
       end do
       
       ! finalize stair
       do i=2, front%n+1
          front%stair(i) = front%stair(i)+front%stair(i-1)
       end do


       ! build rowmap for the rows from the original matrix
       do i=1, front%anrows
          f = first(i) ! f is the front-local column index of the first
          ! coefficient in this row 
          front%stair(f) = front%stair(f)+1
          row = front%stair(f)
          front%arowmap(i) = row
          front%rows(row) = adata%rperm(roff+i)
       end do
       
       ! build rowmap for the rows from the children
       do p = adata%childptr(fnum), adata%childptr(fnum+1)-1
          c = adata%child(p)
          cfront => fdata%front_list(c)
          call qrm_alloc(cfront%rowmap, cfront%ne-cfront%npiv, err)
          __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)
          ! allocate(cfront%rowmap(cfront%ne-cfront%npiv))
          
          ! ne is the number of Householder vectors computed on the
          ! child c. npiv is the number of fully assembled pivots in c
          ne   = cfront%ne
          npiv = cfront%npiv
          
          ! this is the row mapping on the child. cfront%rowmap(k)=i
          ! means that the k-th row of cfront will be assembled into the
          ! i-th row of front
          do i=npiv+1, ne
             f = cfront%colmap(i-npiv)
             front%stair(f) = front%stair(f)+1
             row = front%stair(f)
             cfront%rowmap(i-npiv) = row
          end do
       end do
    end if
    
    ! prevent stair from being above the diagonal. CBs become a nightmare otherwise
    do j=1, front%n
       front%stair(j) = max(front%stair(j),min(j,front%m))
    end do

9999 continue
    if(present(info)) info = err
    deallocate(first)
    return
  end subroutine _qrm_build_stair_rowmapping

  
end subroutine _qrm_activate_front
