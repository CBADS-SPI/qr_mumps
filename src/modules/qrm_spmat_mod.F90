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
!> @file qrm_spmat_mod.F90
!! This file contains the module that implements the main qr_mumps data structure
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


#include "qrm_common.h"


!> @brief This module contains the definition of the basic sparse matrix type
!! and of the associated methods
module _qrm_spmat_mod
  use qrm_common_mod
  use qrm_error_mod
  use qrm_mem_mod
  use qrm_adata_mod
  use _qrm_fdata_mod
  use qrm_error_mod

  !> @brief Generif interface for the @link ::_qrm_spmat_alloc @endlink routine
  interface qrm_spmat_alloc
     module procedure _qrm_spmat_alloc
  end interface

  !> @brief Generif interface for the @link ::_qrm_spmat_init @endlink routine
  interface qrm_spmat_init
     module procedure _qrm_spmat_init
  end interface

  !> @brief Generif interface for the @link ::_qrm_cntl_init @endlink routine
  interface qrm_cntl_init
     module procedure _qrm_cntl_init
  end interface

  !> @brief Generif interface for the @link ::_qrm_spmat_convert @endlink routine
  interface qrm_spmat_convert
     module procedure _qrm_spmat_convert
  end interface

  !> @brief Generif interface for the @link ::_qrm_spmat_copy @endlink routine
  interface qrm_spmat_copy
     module procedure _qrm_spmat_copy
  end interface

  !> @brief Generif interface for the @link ::_qrm_spmat_destroy @endlink routine
  interface qrm_spmat_destroy
     module procedure _qrm_spmat_destroy
  end interface

  !> @brief Generif interface for the @link ::_qrm_spmat_destroy @endlink routine
  interface qrm_spmat_cleanup
     module procedure _qrm_spmat_cleanup
  end interface


  !> @brief Generif interface for the
  !! @link ::_qrm_pseti @endlink,
  !! @link ::_qrm_psetr @endlink and
  ! ! !! @link ::_qrm_psetl @endlink routines
  !!
  interface qrm_set
     module procedure _qrm_pseti, _qrm_psetr!, _qrm_psetl
  end interface

  !> @brief Generif interface for the
  !! @link ::_qrm_pgeti @endlink,
  !! @link ::_qrm_pgetr @endlink and
  !!
  interface qrm_get
     module procedure _qrm_pgeti, _qrm_pgetii, _qrm_pgetr
  end interface

  interface qrm_sync_data
     module procedure :: _qrm_sync_spmat
  end interface qrm_sync_data
  
  interface qrm_sync_spmat
     module procedure  _qrm_sync_spmat
  end interface qrm_sync_spmat

  !> @brief This type defines the data structure used to store a matrix. 

  !> A matrix can be represented either in COO, CSR or CSC
  !! format. Following the qr_mumps convention any array
  !! visible/usable from the users interface will be a pointer;
  !! otherwise allocatbles will be used because they normally provide
  !! better performance.
  !!
  type _qrm_spmat_type
     !> an integer array containing control parameters
     !! expressed as an integer value. Undocumented entries are
     !! assigned as nicntl, nicntl-1 etc.
     !! 
     !! Meaning:
     !! - icntl(qrm_ordering_)     : the ordering method to be used:
     !!             - qrm_auto_    : automatic choice
     !!             - qrm_natural_ : natural order
     !!             - qrm_given_   : given order
     !!             - qrm_colamd_  : COLAMD
     !!             - qrm_metis_   : METIS
     !!             - qrm_scotch_  : SCOTCH
     !! - icntl(qrm_minamalg_)     : minimum node size for amalgamation. Node i can be amalgamated
     !!               to its father iff they both have less than icntl(qrm_minamalg_) variables
     !!               and the fill in introduced is less than rcntl(1).
     !! - icntl(qrm_nb_)           : the block size to be used during factorization
     !! - inctl(qrm_keeph_)        : whether H (the Householder vectors) has to be stored or not 
     !!             - -1  : both R and H are discarded (hidden feature, not documented)
     !!             -  0  : (qrm_no_) only H is discarded
     !!             -  1  : (qrm_yes_) nothing is discarded
     !! - inctl(qrm_ib_)           : the internal blocking size to reduce flops
     !! - icntl(qrm_rhsnb_)        : the blocking parameter to handle multiple right-hand-sides
     !! - icntl(qrm_nlz_)          : the minimum number of subtrees in L0 is this times the
     !!                              number of threads
     !! - icntl(qrm_cnode_)        : the number of cores per NUMA node (or cache, whatever)
     !! - icntl(qrm_bh_)           : in CA trees for the QR facto of a front, this is the number
     !!                              of leaves in each flat subtree
     integer                            :: icntl(nicntl)=0
     !> an double precision array containing control parameters
     !! expressed as a real value. Undocumented entries are
     !! assigned as nrcntl, nrcntl-1 etc.
     !!
     !! Meaning:
     !! - rcntl(qrm_amalgth_=1)      : fill-in threshold for amalgamation. Node i can be amalgamated
     !!                                to its father iff they both have less than icntl(3) variables
     !!                                and the fill in introduced is less than rcntl(1).
     !! - rcntl(qrm_rweight_=nrcntl) : subtrees with relative weight below this threashold will
     !!                                be put in L0
     !! - rcntl(qrm_mem_relax_=2)    : the relaxation parameter for the factorization memory
     real(kind(1.d0))                   :: rcntl(nrcntl)=0.d0
     !> an array containing local global stats. some of its content will only
     !! be relevant on the master node. Meaning:
     !! - gstats(1) = estimated total number of flops
     !! - gstats(2) = estimated total number of nonzeroes in R (estimated)
     !! - gstats(3) = estimated total number of nonzeroes in H (estimated)
     !! - gstats(4) = actual total number of flops
     !! - gstats(5) = actual total number of nonzeroes in R (estimated)
     !! - gstats(6) = actual total number of nonzeroes in H (estimated)
     !! - gstats(7) = estimated peak memory consumption
     integer(kind=8)                    :: gstats(10)=0
     !> Pointer to the beginning of each row in CSR format
     integer, pointer,    dimension(:)  :: iptr => null()
     !> Pointer to the beginning of each column in CSC format 
     integer, pointer,    dimension(:)  :: jptr => null()
     !> Row indices
     integer, pointer,    dimension(:)  :: irn => null()
     !> Column indices
     integer, pointer,    dimension(:)  :: jcn => null()
     !> Numerical values
     _qrm_data, pointer, dimension(:)   :: val => null()
     !> Number of rows
     integer                            :: m=0
     !> Number of columns
     integer                            :: n=0
     !> Number of nonzero elements
     integer                            :: nz=0
     !> A pointer to an array containing a column permutation provided
     !! by the user
     integer, pointer,    dimension(:)  :: cperm_in => null()
     !> a @link qrm_adata_mod::qrm_adata_type @endlink data which is meant to
     !> contain all the data related to the analysis phase         
     type(qrm_adata_type), allocatable  :: adata
     !> a @link _qrm_fdata_mod::_qrm_fdata_type @endlink data which is meant to
     !> contain all the data related to the factorization phase         
     type(_qrm_fdata_type), allocatable :: fdata
     !> Storage format; can be either 'COO', 'CSR' or 'CSC'
     character(len=3)                   :: fmt='coo'
  end type _qrm_spmat_type

  
contains

  !> @brief This subroutine allocates memory for a sparse matrix.
  !!
  !! @param[in,out] qrm_spmat A @link _qrm_spmat_mod::_qrm_spmat_type @endlink data
  !!           structure. The memory for storing the matrix is
  !!           allocated according to the storage format. Also
  !!           qrm_spmat%nz, qrm_spmat%m and qrm_spmat%n are set to
  !!           nz, m and n respectively.
  !!           These are the sizes of the arrays in output
  !!           * coo: irn(nz), jcn(nz), val(nz)
  !!           * csr: iptr(m+1), jcn(nz), val(nz)
  !!           * csc: irn(nz), jptr(n+1), val(nz)
  !!
  !! @param[in] nz  The number of nonzeroes contained in the matrix
  !!
  !! @param[in] m   The number of rows in the matrix
  !!
  !! @param[in] n   The number of columns in the matrix
  !!
  !! @param[in] fmt The matrix storage format. Can be either "coo" or "csr"
  !!           or "csc"
  !!
  subroutine _qrm_spmat_alloc(qrm_spmat, nz, m, n, fmt, info)
    use qrm_mem_mod
    implicit none

    type(_qrm_spmat_type), intent(inout) :: qrm_spmat
    integer, intent(in)                  :: nz, m, n
    character, intent(in)                :: fmt*(*)
    integer, optional                    :: info

    integer                              :: err
    character(len=*), parameter          :: name='qrm_spmat_alloc'

    err = 0
    
#if defined(debug)
    __QRM_PRNT_DBG('("Allocating Matrix")')
#endif

    if(fmt .eq. 'coo') then
       if(err.eq.0) call qrm_alloc(qrm_spmat%irn, nz, info=err)
       if(err.eq.0) call qrm_alloc(qrm_spmat%jcn, nz, info=err)
       if(err.eq.0) call qrm_alloc(qrm_spmat%val, nz, info=err)
       __QRM_INFO_CHECK(err, name, 'qrm_alloc', 9999)
    else if(fmt .eq. 'csr') then
       if(err.eq.0) call qrm_alloc(qrm_spmat%iptr, m+1, info=err)
       if(err.eq.0) call qrm_alloc(qrm_spmat%jcn , nz,  info=err)
       if(err.eq.0) call qrm_alloc(qrm_spmat%val , nz,  info=err)
       __QRM_INFO_CHECK(err, name, 'qrm_alloc', 9999)
    else if(fmt .eq. 'csc') then
       if(err.eq.0) call qrm_alloc(qrm_spmat%irn , nz,  info=err)
       if(err.eq.0) call qrm_alloc(qrm_spmat%jptr, n+1, info=err)
       if(err.eq.0) call qrm_alloc(qrm_spmat%val , nz,  info=err)
       __QRM_INFO_CHECK(err, name, 'qrm_alloc', 9999)
    else
       ! format not known. set info and jump to 9999
       err = 1
       call qrm_error_print(err, name, aed=fmt)
       goto 9999
    end if

    qrm_spmat%nz = nz
    qrm_spmat%m  = m
    qrm_spmat%n  = n

9999 continue
    if(present(info)) info = err
    return

  end subroutine _qrm_spmat_alloc


  !> @brief This subroutine initializes a qrm_spmat_type instance setting
  !! default values into the control parameters
  !!
  !! @param[in,out] qrm_spmat The matrix to be initialized
  !! 
  subroutine _qrm_spmat_init(qrm_spmat, info)

    implicit none
    type(_qrm_spmat_type), intent(inout) :: qrm_spmat
    integer, optional                    :: info

    
    call _qrm_cntl_init(qrm_spmat)

    nullify(qrm_spmat%iptr, qrm_spmat%jptr, qrm_spmat%irn, qrm_spmat%jcn, &
         & qrm_spmat%val, qrm_spmat%cperm_in)

    if(present(info)) info = 0
    return

  end subroutine _qrm_spmat_init


  !> @brief This subroutine initializes a qrm_spmat_type instance setting
  !! default values into the control parameters
  !!
  !! @param[in,out] qrm_spmat The matrix to be initialized
  !! 
  subroutine _qrm_cntl_init(qrm_spmat)

    use qrm_common_mod
    implicit none
    type(_qrm_spmat_type), intent(inout) :: qrm_spmat

    ! set default values for icntl and rcntl
    qrm_spmat%icntl(qrm_ordering_)      = qrm_auto_
    qrm_spmat%icntl(qrm_minamalg_)      = 4
    qrm_spmat%icntl(qrm_mb_)            = 256
    qrm_spmat%icntl(qrm_nb_)            = 128
    qrm_spmat%icntl(qrm_ib_)            = 32
    qrm_spmat%icntl(qrm_bh_)            = -1
    qrm_spmat%icntl(qrm_keeph_)         = qrm_yes_
    qrm_spmat%icntl(qrm_rhsnb_)         = -1
    qrm_spmat%icntl(qrm_nlz_)           = 8
    qrm_spmat%icntl(qrm_cnode_)         = 1
    qrm_spmat%icntl(qrm_sing_)          = qrm_no_
    
    qrm_spmat%rcntl(qrm_amalgth_)       = 0.05
    qrm_spmat%rcntl(qrm_rweight_)       = 0.001
    qrm_spmat%rcntl(qrm_mem_relax_)     = -100.0
    qrm_spmat%fmt = 'coo'

    return

  end subroutine _qrm_cntl_init





  !> This subroutine converts an input matrix into a different
  !! storage format. Optionally the values may be ignored
  !! (this comes handy during the analysis)
  !!
  !! @param[in] in_mat      the input matrix
  !!
  !! @param[in,out] out_mat the output matrix in fmt format
  !!
  !! @param[in] fmt         the format of the output matrix
  !!
  !! @param[in]  values      (optional) if values=.true. the output matrix will include
  !!                        numerical values, otherwise only the structure
  !!
  subroutine _qrm_spmat_convert(in_mat, out_mat, fmt, values, info)
    implicit none

    type(_qrm_spmat_type), intent(in)  :: in_mat
    type(_qrm_spmat_type)              :: out_mat
    character, intent(in)              :: fmt*(*)
    logical, optional                  :: values
    integer, optional                  :: info

    integer                            :: err
    character(len=*), parameter        :: name='qrm_spmat_convert'

    err = 0
    
    select case(in_mat%fmt)
    case('csc')
       select case(fmt)
       case('csr')
          call _qrm_csc_to_csr(in_mat, out_mat, values, err)
          __QRM_INFO_CHECK(err, name, 'qrm_csc_to_csr', 9999)
       case default
          ! format not known. set info and jump to 9999
          err = 1
          call qrm_error_print(err, name, aed=fmt)
          goto 9999
       end select
    case('coo')
       select case(fmt)
       case('csc')
          call _qrm_coo_to_csc(in_mat, out_mat, values, err)
          __QRM_INFO_CHECK(err, name, 'qrm_coo_to_csc', 9999)
       case('csr')
          call _qrm_coo_to_csr(in_mat, out_mat, values, err)
          __QRM_INFO_CHECK(err, name, 'qrm_coo_to_csr', 9999)
       case default
          ! format not known. set info and jump to 9999
          err = 1
          call qrm_error_print(err, name, aed=fmt)
          goto 9999
       end select
    case default
       ! format not known. set info and jump to 9999
       info = 1
       call qrm_error_print(info, name, aed=fmt)
       goto 9999
    end select
    
    out_mat%icntl    = in_mat%icntl
    out_mat%rcntl    = in_mat%rcntl

9999 continue
    if(present(info)) info = err
    return

  end subroutine _qrm_spmat_convert


  !> This subroutine converts a COO matrix into a CSC
  !! matrix. Optionally the values may be ignored
  !! (this comes handy during the analysis)
  !!
  !! @param[in] in_mat      the input matrix
  !!
  !! @param[in,out] out_mat the output matrix in fmt format
  !!
  !! @param[in]  values      (optional) if values=.true. the output matrix will include
  !!                        numerical values, otherwise only the structure
  !!
  subroutine _qrm_coo_to_csc(in_mat, out_mat, values, info)

    implicit none

    type(_qrm_spmat_type), intent(in)  :: in_mat
    type(_qrm_spmat_type)              :: out_mat
    logical, optional                  :: values
    integer, optional                  :: info

    integer, allocatable               :: work(:)
    logical                            :: ivalues, ob
    integer                            :: i, j, idx, k, m, n

    ! error management
    integer                            :: err
    character(len=*), parameter        :: name='_qrm_coo_to_csc'

    err = 0
    
    if(present(values)) then
       ivalues = values
    else
       ivalues = .true.
    end if

    if(ivalues)  call qrm_realloc(out_mat%val,  in_mat%nz,  err)
    if(err.eq.0) call qrm_alloc(work,           in_mat%n+1, err)
    if(err.eq.0) call qrm_realloc(out_mat%jptr, in_mat%n+1, err)
    if(err.eq.0) call qrm_realloc(out_mat%irn,  in_mat%nz,  err)
    __QRM_INFO_CHECK(err, name, 'qrm_(re)alloc', 9999)

    work=0
    ob = .false.

    m = in_mat%m
    n = in_mat%n

    ! first loop to calculate # of nz per column
    do k=1, in_mat%nz
       j = in_mat%jcn(k)
       i = in_mat%irn(k)
       if((j.gt.0) .and. (j.le. n) .and. (i.gt.0) .and. (i.le. m) ) then
          work(j) = work(j)+1
       else
          ! out of bounds coefficients. ignore and print a warning at the end
          ob = .true.
       end if
    end do

    if(ob) then
       __QRM_PRNT_DBG('("** Out-of-bounds coefficients present **")')
    end if
    
    ! loop to convert the counts into ptrs
    out_mat%jptr(1) = 1
    do j=2, n+1
       out_mat%jptr(j) = out_mat%jptr(j-1)+work(j-1)
    end do


    ! last loop to put things into place
    work=0
    ! instead of putting an "if" inside the loop
    ! I put it here to gain some speed
    if(ivalues) then
       do k=1, in_mat%nz
          j = in_mat%jcn(k)
          i = in_mat%irn(k)
          if((j.le.0) .or. (j.gt. n) .or. (i.le.0) .or. (i.gt. m) ) cycle
          idx = out_mat%jptr(j)+work(j)
          out_mat%irn(idx) = i
          out_mat%val(idx) = in_mat%val(k)
          work(j) = work(j)+1
       end do
    else
       do k=1, in_mat%nz
          j = in_mat%jcn(k)
          i = in_mat%irn(k)
          if((j.le.0) .or. (j.gt. n) .or. (i.le.0) .or. (i.gt. m) ) cycle
          idx = out_mat%jptr(j)+work(j)
          out_mat%irn(idx) = i
          work(j) = work(j)+1
       end do
    end if

    out_mat%m   = in_mat%m
    out_mat%n   = in_mat%n
    out_mat%nz  = in_mat%nz
    out_mat%fmt = 'csc'
    
9999 continue
    ! cleanup and return
    call qrm_dealloc(work, info)

    if(err.ne.0) then
       call qrm_dealloc(out_mat%jptr)
       call qrm_dealloc(out_mat%irn)
       call qrm_dealloc(out_mat%val)
    end if

    if(present(info)) info = err
    return

  end subroutine _qrm_coo_to_csc


  !> This subroutine converts a COO matrix into a CSC
  !! matrix. Optionally the values may be ignored
  !! (this comes handy during the analysis)
  !!
  !! @param[in] in_mat      the input matrix
  !!
  !! @param[in,out] out_mat the output matrix in fmt format
  !!
  !! @param[in]  values      (optional) if values=.true. the output matrix will include
  !!                        numerical values, otherwise only the structure
  !!
  subroutine _qrm_coo_to_csr(in_mat, out_mat, values, info)

    implicit none

    type(_qrm_spmat_type), intent(in)  :: in_mat
    type(_qrm_spmat_type)              :: out_mat
    logical, optional                  :: values
    integer, optional                  :: info

    integer, allocatable               :: work(:)
    logical                            :: ivalues, ob
    integer                            :: i, j, idx, k, m, n

    ! error management
    integer                            :: err
    character(len=*), parameter        :: name='_qrm_coo_to_csr'

    err = 0

    if(present(values)) then
       ivalues = values
    else
       ivalues = .true.
    end if

    if(ivalues) call qrm_realloc(out_mat%val , in_mat%nz, err)
    if(err.eq.0) call qrm_alloc(work, in_mat%m+1, err)
    if(err.eq.0) call qrm_realloc(out_mat%iptr, in_mat%m+1, err)
    if(err.eq.0) call qrm_realloc(out_mat%jcn , in_mat%nz, err)
    __QRM_INFO_CHECK(err, name,'qrm_(re)alloc',9999)

    work=0
    ob = .false.

    m = in_mat%m
    n = in_mat%n

    ! first loop to calculate # of nz per row
    do k=1, in_mat%nz
       j = in_mat%jcn(k)
       i = in_mat%irn(k)
       if((i.gt.0) .and. (i.le. m) .and. (j.gt.0) .and. (j.le.n) ) then
          work(i) = work(i)+1
       else
          ! out of bounds coefficients. ignore and print a warning at the end
          ob = .true.
       end if
    end do

    if(ob) then
       __QRM_PRNT_DBG('("** Out-of-bounds coefficients present **")')
    end if
    
    ! loop to convert the counts into ptrs
    out_mat%iptr(1) = 1
    do i=2, m+1
       out_mat%iptr(i) = out_mat%iptr(i-1)+work(i-1)
    end do


    ! last loop to put things into place
    work=0
    ! instead of putting an "if" inside the loop
    ! I put it here to gain some speed
    if(ivalues) then
       do k=1, in_mat%nz
          j = in_mat%jcn(k)
          i = in_mat%irn(k)
          if((j.le.0) .or. (j.gt. n) .or. (i.le.0) .or. (i.gt. m) ) cycle
          idx = out_mat%iptr(i)+work(i)
          out_mat%jcn(idx) = j
          out_mat%val(idx) = in_mat%val(k)
          work(i) = work(i)+1
       end do
    else
       do k=1, in_mat%nz
          j = in_mat%jcn(k)
          i = in_mat%irn(k)
          if((j.le.0) .or. (j.gt. n) .or. (i.le.0) .or. (i.gt. m) ) cycle
          idx = out_mat%iptr(i)+work(i)
          out_mat%jcn(idx) = j
          work(i) = work(i)+1
       end do
    end if

    out_mat%m   = in_mat%m
    out_mat%n   = in_mat%n
    out_mat%nz  = in_mat%nz
    out_mat%fmt = 'csr'

9999 continue 
    ! cleanup and return
    call qrm_dealloc(work)

    if(err.ne.0) then
       call qrm_dealloc(out_mat%iptr)
       call qrm_dealloc(out_mat%jcn)
       call qrm_dealloc(out_mat%val)
    end if
    
    if(present(info)) info = err
    return

  end subroutine _qrm_coo_to_csr

  
  !> This subroutine converts a CSC matrix into a CSR
  !! matrix. Optionally the values may be ignored
  !! (this comes handy during the analysis)
  !!
  !! @param[in] in_mat      the input matrix
  !!
  !! @param[in,out] out_mat the output matrix in fmt format
  !!
  !! @param[in]  values      (optional) if values=.true. the output matrix will include
  !!                        numerical values, otherwise only the structure
  !!
  subroutine _qrm_csc_to_csr(in_mat, out_mat, values, info)
    implicit none
    
    type(_qrm_spmat_type), intent(in)  :: in_mat
    type(_qrm_spmat_type)              :: out_mat
    logical, optional                  :: values
    integer, optional                  :: info

    integer, allocatable               :: work(:)

    logical                            :: ivalues, ob
    integer                            :: i, j, idx, ii, m, n

    ! error management
    integer                            :: err
    character(len=*), parameter        :: name='_qrm_csc_to_csr'

    err = 0

    if(present(values)) then
       ivalues=values
    else
       ivalues = .true.
    end if

    ob = .false.

    m = in_mat%m
    n = in_mat%n

    if(ivalues)  call qrm_realloc(out_mat%val,  in_mat%nz, err)
    if(err.eq.0) call qrm_alloc(work,           m+1,       err)
    if(err.eq.0) call qrm_realloc(out_mat%iptr, m+1,       err)
    if(err.eq.0) call qrm_realloc(out_mat%jcn,  in_mat%nz, err)
    __QRM_INFO_CHECK(err, name,'qrm_(re)alloc',9999)

    work=0
    ! first loop to calculate # of nz per row
    do j = 1, n
       do ii= in_mat%jptr(j), in_mat%jptr(j+1)-1
          i = in_mat%irn(ii)
          if((i.gt.0) .and. (i.le.m)) then
             work(i) = work(i)+1
          else
             ob = .true.
          end if
       end do
    end do

    if(ob) then
       __QRM_PRNT_DBG('("** Out-of-bounds coefficients present **")')
    end if

    ! loop to convert the counts into ptrs
    out_mat%iptr(1) = 1
    do j=2, m+1
       out_mat%iptr(j) = out_mat%iptr(j-1)+work(j-1)
    end do


    ! last loop to put things into place
    work=0
    ! instead of putting an "if" inside the loop
    ! I put it here to gain some speed
    if(ivalues) then
       do j = 1, n
          do ii= in_mat%jptr(j), in_mat%jptr(j+1)-1
             i = in_mat%irn(ii)
             if((i.le.0) .or. (i.gt.m)) cycle
             idx = out_mat%iptr(i)+work(i)
             out_mat%jcn(idx) = j
             out_mat%val(idx) = in_mat%val(ii)
             work(i) = work(i)+1
          end do
       end do
    else
       do j = 1, n
          do ii= in_mat%jptr(j), in_mat%jptr(j+1)-1
             i = in_mat%irn(ii)
             if((i.le.0) .or. (i.gt.m)) cycle
             idx = out_mat%iptr(i)+work(i)
             out_mat%jcn(idx) = j
             work(i) = work(i)+1
          end do
       end do
    end if

    out_mat%m   = in_mat%m
    out_mat%n   = in_mat%n
    out_mat%nz  = in_mat%nz
    out_mat%fmt = 'csr'

9999 continue
    ! cleanup and return
    call qrm_dealloc(work)
    if(err.ne.0) then
       call qrm_dealloc(out_mat%iptr)
       call qrm_dealloc(out_mat%jcn)
       call qrm_dealloc(out_mat%val)
    end if
    
    if(present(info)) info = err
    return

  end subroutine _qrm_csc_to_csr


  !> This subroutine sorts the values of a matrix in order of
  !> increasing indexes
  !!
  !! @param[in] qrm_mat      the input matrix to be sorted
  !!
  subroutine _qrm_sort_mat(qrm_mat, values, rc, info)

    implicit none
    type(_qrm_spmat_type)       :: qrm_mat
    logical, optional           :: values
    character(len=3), optional  :: rc
    integer, optional           :: info

    character(len=3)            :: irc
    logical                     :: ivalues

    ! error management
    integer                     :: err
    character(len=*), parameter :: name='_qrm_sort_mat'
    
    err = 0

    if(present(rc)) then
       irc = rc
    else
       irc = 'row'
    end if

    if(present(values)) then
       ivalues = values
    else
       ivalues = .true.
    end if

    select case(qrm_mat%fmt)
    case('csc')
       call _qrm_sort_csc_mat(qrm_mat, ivalues, err)
       __QRM_INFO_CHECK(err, name,'qrm_sort_csc_mat',9999)
    case('csr')
       call _qrm_sort_csr_mat(qrm_mat, ivalues, err)
       __QRM_INFO_CHECK(err, name,'qrm_sort_csr_mat',9999)
    case('coo')
       call _qrm_sort_coo_mat(qrm_mat, ivalues, irc, err)
       __QRM_INFO_CHECK(err, name,'qrm_sort_coo_mat',9999)
    end select

9999 continue
    if(present(info)) info = err
    return
  end subroutine _qrm_sort_mat


  !> This subroutine sorts the values of a CSC matrix in order of
  !> increasing indexes
  !!
  !! @param[in] qrm_mat      the input matrix to be sorted
  !!
  subroutine _qrm_sort_csc_mat(qrm_mat, values, info)
    use qrm_sort_mod
    implicit none
    type(_qrm_spmat_type)       :: qrm_mat
    logical                     :: values
    integer, optional           :: info

    integer, allocatable        :: aux(:)
    integer                     :: j, m, k1, k2

    integer                     :: err
    character(len=*), parameter :: name='_qrm_sort_csc_mat'
    
    err = 0

    call qrm_alloc(aux, qrm_mat%m+2, err)
    __QRM_INFO_CHECK(err, name, 'qrm_alloc',9999)
    
    do j=1,qrm_mat%m
       k1= qrm_mat%jptr(j)
       k2= qrm_mat%jptr(j+1)-1
       m = k2-k1+1
       call qrm_mergesort(m, qrm_mat%irn(k1:k2), aux(1:m+2))
       if(values) then
          call qrm_mergeswap(m, aux(1:m+2), qrm_mat%irn(k1:k2), qrm_mat%val(k1:k2))
       else
          call qrm_mergeswap(m, aux(1:m+2), qrm_mat%irn(k1:k2))
       end if
    end do

    call qrm_dealloc(aux, info)

9999 continue
    if(present(info)) info = err
    return
  end subroutine _qrm_sort_csc_mat


  !> This subroutine sorts the values of a CSR matrix in order of
  !> increasing indexes
  !!
  !! @param[in] qrm_mat      the input matrix to be sorted
  !!
  subroutine _qrm_sort_csr_mat(qrm_mat, values, info)
    use qrm_sort_mod
    implicit none

    type(_qrm_spmat_type)       :: qrm_mat
    logical                     :: values
    integer, optional           :: info

    integer, allocatable        :: aux(:)
    integer                     :: i, n, k1, k2

    integer                     :: err
    character(len=*), parameter :: name='_qrm_sort_csr_mat'

    err = 0

    call qrm_alloc(aux, qrm_mat%n+2, err)
    __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)

    do i=1,qrm_mat%m
       k1= qrm_mat%iptr(i)
       k2= qrm_mat%iptr(i+1)-1
       n = k2-k1+1
       call qrm_mergesort(n, qrm_mat%jcn(k1:k2), aux(1:n+2))
       if(values) then
          call qrm_mergeswap(n, aux(1:n+2), qrm_mat%jcn(k1:k2), qrm_mat%val(k1:k2))
       else
          call qrm_mergeswap(n, aux(1:n+2), qrm_mat%jcn(k1:k2))
       end if
    end do

    call qrm_dealloc(aux)

9999 continue

    if(present(info)) info = err
    return
    
  end subroutine _qrm_sort_csr_mat


  !> This subroutine sorts the values of a COO matrix in order of
  !> increasing indexes
  !!
  !! @param[in] qrm_mat      the input matrix to be sorted
  !!
  subroutine _qrm_sort_coo_mat(qrm_mat, values, rc, info)
    use qrm_sort_mod
    implicit none
    type(_qrm_spmat_type)       :: qrm_mat
    logical                     :: values
    character(len=3)            :: rc
    integer, optional           :: info

    integer, allocatable        :: aux(:)
    integer                     :: j, m, k1, k2, nz

    integer                     :: err
    character(len=*), parameter :: name="qrm_sort_coo_mat"
    
    err = 0

    nz = qrm_mat%nz
    
    call qrm_alloc(aux, nz+2, err)
    __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)

    if(rc.eq.'row') then
       call qrm_mergesort(nz, qrm_mat%irn(1:nz), aux(1:nz+2))
    else if(rc.eq.'col') then
       call qrm_mergesort(nz, qrm_mat%jcn(1:nz), aux(1:nz+2))
    end if

    if(values) then
       call qrm_mergeswap(nz, aux(1:nz+2), qrm_mat%irn(1:nz), qrm_mat%jcn(1:nz), qrm_mat%val(1:nz))
    else
       call qrm_mergeswap(nz, aux(1:nz+2), qrm_mat%irn(1:nz), qrm_mat%jcn(1:nz))
    end if

    k1 = 1
    if(rc.eq.'row') then
       outerr: do
          if (k1.ge.qrm_mat%nz) exit outerr
          k2 = k1
          innerr: do
             if(k2.eq.qrm_mat%nz) exit innerr
             if(qrm_mat%irn(k2+1).ne.qrm_mat%irn(k1)) exit
             k2 = k2+1
          end do innerr
          m = k2-k1+1
          call qrm_mergesort(m, qrm_mat%jcn(k1:k2), aux(1:m+2))
          if(values) then
             call qrm_mergeswap(m, aux(1:m+2), qrm_mat%jcn(k1:k2), qrm_mat%val(k1:k2))
          else
             call qrm_mergeswap(m, aux(1:m+2), qrm_mat%jcn(k1:k2))
          end if
          k1 = k2+1
       end do outerr
    else if(rc.eq.'col') then
       outerc: do
          if (k1.ge.qrm_mat%nz) exit outerc
          k2 = k1
          innerc: do
             if(k2.eq.qrm_mat%nz) exit innerc
             if(qrm_mat%jcn(k2+1).ne.qrm_mat%jcn(k1)) exit
             k2 = k2+1
          end do innerc
          m = k2-k1+1
          call qrm_mergesort(m, qrm_mat%irn(k1:k2), aux(1:m+2))
          if(values) then
             call qrm_mergeswap(m, aux(1:m+2), qrm_mat%irn(k1:k2), qrm_mat%val(k1:k2))
          else
             call qrm_mergeswap(m, aux(1:m+2), qrm_mat%irn(k1:k2))
          end if
          k1 = k2+1
       end do outerc
    end if

    call qrm_dealloc(aux)

9999 continue

    if(present(info)) info = err
    return
  end subroutine _qrm_sort_coo_mat
  

  !> This subroutine makes a copy of a matrix. Optionally the values
  !! may be ignored (this comes handy during the analysis)
  !!
  !! @param[in] in_mat      the input matrix
  !!
  !! @param[in,out] out_mat the output matrix in fmt format
  !!
  !! @param[in]  values      (optional) if values=.true. the output matrix will include
  !!                        numerical values, otherwise only the structure
  !!
  subroutine _qrm_spmat_copy(in_mat, out_mat, values, info)

    type(_qrm_spmat_type), intent(in) :: in_mat
    type(_qrm_spmat_type)             :: out_mat
    logical, optional                 :: values
    integer, optional                 :: info

    logical                           :: ivalues=.true.
    ! error management
    integer                           :: err
    character(len=*), parameter       :: name='_qrm_spmat_copy'

    err = 0

    if(present(values)) then
       ivalues=values
    else
       ivalues=.true.
    end if

    select case(in_mat%fmt)
    case('csc')
       if(ivalues)  call qrm_realloc(out_mat%val,  in_mat%nz,  err)
       if(err.eq.0) call qrm_realloc(out_mat%jptr, in_mat%n+1, err)
       if(err.eq.0) call qrm_realloc(out_mat%irn,  in_mat%nz,  err)
       __QRM_INFO_CHECK(err, name,'qrm_realloc',9999)

       if(ivalues) out_mat%val(1:in_mat%nz) = in_mat%val(1:in_mat%nz)
       out_mat%jptr(1:in_mat%n+1)           = in_mat%jptr(1:in_mat%n+1)
       out_mat%irn(1:in_mat%nz)             = in_mat%irn(1:in_mat%nz)
    case('coo')
       if(ivalues)  call qrm_realloc(out_mat%val, in_mat%nz, err)
       if(err.eq.0) call qrm_realloc(out_mat%jcn, in_mat%nz, err)
       if(err.eq.0) call qrm_realloc(out_mat%irn, in_mat%nz, err)
       __QRM_INFO_CHECK(err, name,'qrm_realloc',9999)

       if(ivalues) out_mat%val(1:in_mat%nz) = in_mat%val(1:in_mat%nz)
       out_mat%jcn(1:in_mat%nz)             = in_mat%jcn(1:in_mat%nz)
       out_mat%irn(1:in_mat%nz)             = in_mat%irn(1:in_mat%nz)
    case default
       ! format not known. set info and jump to 9999
       info = 1
       call qrm_error_print(info, name, aed=in_mat%fmt)
       goto 9999
    end select

    out_mat%n        = in_mat%n
    out_mat%m        = in_mat%m
    out_mat%nz       = in_mat%nz
    out_mat%fmt      = in_mat%fmt
    out_mat%icntl    = in_mat%icntl
    out_mat%rcntl    = in_mat%rcntl

9999 continue
    ! cleanup and return
    if(err.ne.0) call _qrm_spmat_destroy(out_mat, all=.true.)

    if(present(info)) info = err
    return

  end subroutine _qrm_spmat_copy

  !> @brief This subroutine destroyes a qrm_spmat instance
  !!
  !! @param[in,out] qrm_spmat the matrix to be destroyed
  !!
  !! @param[in] all whether to deallocate all the memory or not
  !!
  subroutine _qrm_spmat_destroy(qrm_spmat, all, info)

    use qrm_mem_mod
    implicit none

    type(_qrm_spmat_type)           :: qrm_spmat
    logical, optional               :: all
    integer, optional               :: info

    logical                         :: iall
    ! error management
    integer                         :: err
    character(len=*), parameter     :: name='_qrm_spmat_destroy'

    err = 0

    if(present(all)) then
       iall = all
    else
       iall = .false.
    end if

    if(iall) then
       call qrm_dealloc(qrm_spmat%irn,      info=err)
       call qrm_dealloc(qrm_spmat%jcn,      info=err)
       call qrm_dealloc(qrm_spmat%iptr,     info=err)
       call qrm_dealloc(qrm_spmat%jptr,     info=err)
       call qrm_dealloc(qrm_spmat%val,      info=err)
       call qrm_dealloc(qrm_spmat%cperm_in, info=err)
       __QRM_INFO_CHECK(err, name, 'qrm_dealloc', 9999)
    end if

    qrm_spmat%n       = 0
    qrm_spmat%m       = 0
    qrm_spmat%nz      = 0
    qrm_spmat%fmt     = ''

    call _qrm_spmat_cleanup(qrm_spmat, err)
    __QRM_INFO_CHECK(err, name, 'qrm_spmat_cleanup', 9999)

9999 continue

    if(present(info)) info = err
    return

  end subroutine _qrm_spmat_destroy


  !> @brief This subroutine destroyes a qrm_spmat instance
  !!
  !! @param[in,out] qrm_spmat the matrix to be destroyed
  !!
  !! @param[in] all whether to deallocate all the memory or not
  !!
  subroutine _qrm_spmat_cleanup(qrm_spmat, info)

    use qrm_mem_mod
    implicit none

    type(_qrm_spmat_type)       :: qrm_spmat
    integer, optional           :: info
    
    ! error management
    integer                     :: err
    character(len=*), parameter :: name='_qrm_spmat_cleanup'

    err = 0
    
    if(allocated(qrm_spmat%adata)) then
       call qrm_adata_destroy(qrm_spmat%adata, err)
       __QRM_INFO_CHECK(err, name, 'qrm_adata_destroy', 9999)
    end if
    if(allocated(qrm_spmat%fdata)) then
       call _qrm_fdata_destroy(qrm_spmat%fdata, err)
       __QRM_INFO_CHECK(err, name, 'qrm_fdata_destroy', 9999)
    end if

9999 continue

    if(present(info)) info = err
    return

  end subroutine _qrm_spmat_cleanup



  ! The following subroutine set or get control parameters from the
  ! cntl or rcntl control arrays. All the set and get routines are
  ! gathered under the same, overloaded interface, respectively


  !> @brief This subroutine is meant to set the integer control parameters
  !!
  !! @param[in,out] qrm_spmat The qrm_spmat instance concerned by the setting
  !!
  !! @param[in] string a string describing the parameter to be
  !!            set. Accepted values are:
  !!            - "qrm_ordering" : to set a method for the fill-reducing
  !!                               column permutation. Accepted values are:
  !!                               - qrm_auto_=0    for automatic choice
  !!                               - qrm_natural_=1 for natural ordering
  !!                               - qrm_given_=2   for given ordering (through
  !!                                                the qrm_spmat%cperm_in pointer)
  !!                               - qrm_colamd_=3  for COLAMD
  !!                               - qrm_metis_=4   for METIS
  !!                               - qrm_scotch_=5  for SCOTCH
  !!            - "qrm_minamalg" : fronts whose size is smaller than this will be
  !!                               systematically amalgamated to their father
  !!            - "qrm_nb" : the block-size that defines the granularity of parallel tasks
  !!            - "qrm_ib" : the block-size for computations
  !!            - "qrm_rhsnb" : the block-size for grouping RHSs
  !!            - "qrm_keeph" : whether to store or not the Householder vectors. Accepted
  !!                            values are qrm_yes_ and qrm_no_
  !!            - "qrm_sing" : whether or not to detect the presence of singletons. Accepted
  !!                           values are qrm_yes_ and qrm_no_
  !!            - "qrm_nlz"  : the number of subtrees in L0 will be at least this times
  !!                           the number of threads
  !!            - "qrm_cnode" : the number of cores per node
  !!
  !! @param[in] ival Any of the accepted values described above
  !!
  subroutine _qrm_pseti(qrm_spmat, string, ival, info)
    use qrm_common_mod
    use qrm_string_mod
    use qrm_error_mod
    implicit none

    type(_qrm_spmat_type)           :: qrm_spmat
    character(len=*)                :: string
    integer                         :: ival
    integer, optional               :: info

    character(len=len(string))      :: istring

    ! error management
    integer                         :: err
    character(len=*), parameter     :: name='_qrm_pseti'

    err = 0
    
    istring = qrm_str_tolower(string)
    if(index(istring,'qrm_ordering') .eq. 1) then
       qrm_spmat%icntl(qrm_ordering_) = ival
    else if (index(istring,'qrm_minamalg') .eq. 1) then
       qrm_spmat%icntl(qrm_minamalg_) = ival
       qrm_spmat%icntl(qrm_nb_) = ival
    else if (index(istring,'qrm_mb') .eq. 1) then
       qrm_spmat%icntl(qrm_mb_) = ival
    else if (index(istring,'qrm_nb') .eq. 1) then
       qrm_spmat%icntl(qrm_nb_) = ival
    else if (index(istring,'qrm_ib') .eq. 1) then
       qrm_spmat%icntl(qrm_ib_) = ival
    else if (index(istring,'qrm_bh') .eq. 1) then
       qrm_spmat%icntl(qrm_bh_) = ival
    else if (index(istring,'qrm_rhsnb') .eq. 1) then
       qrm_spmat%icntl(qrm_rhsnb_) = ival
    else if (index(istring,'qrm_keeph') .eq. 1) then
       select case (ival)
       case (:-1)
          qrm_spmat%icntl(qrm_keeph_) = -1
       case (0)
          qrm_spmat%icntl(qrm_keeph_) = 0
       case (1:)
          qrm_spmat%icntl(qrm_keeph_) = 1
       end select
    else if (index(istring,'qrm_sing') .eq. 1) then
       if(ival .eq. qrm_yes_) then
          qrm_spmat%icntl(qrm_sing_) = ival
       else
          qrm_spmat%icntl(qrm_sing_) = qrm_no_
       end if
    else if (index(istring,'qrm_nlz') .eq. 1) then
       qrm_spmat%icntl(qrm_nlz_) = ival
    else if (index(istring,'qrm_cnode') .eq. 1) then
       qrm_spmat%icntl(qrm_cnode_) = ival
    else
       err = 23
       call qrm_error_print(err, name, aed=string)
       goto 9999
    end if

9999 continue
    if(present(info)) info = err
    return

  end subroutine _qrm_pseti


  !> @brief This subroutine is meant to set the real control parameters
  !!
  !! @param[in,out] qrm_spmat The qrm_spmat instance concerned by the setting
  !!
  !! @param[in] string a string describing the parameter to be
  !!            set. Accepted values are:
  !!            - "qrm_amalgth" : the threshold that controls the amalgamation.
  !!                              A higher threshold means more fill-in but also
  !!                              more BLAS-3. Any real value is accepted
  !!
  !! @param[in] rval Any of the accepted values described above
  !!
  subroutine _qrm_psetr(qrm_spmat, string, rval, info)

    use qrm_common_mod
    use qrm_string_mod
    use qrm_error_mod
    implicit none

    type(_qrm_spmat_type)           :: qrm_spmat
    character(len=*)                :: string
    real(kind(1.d0))                :: rval
    integer, optional               :: info

    character(len=len(string))      :: istring

    ! error management
    integer                         :: err
    character(len=*), parameter     :: name='_qrm_psetr'

    err = 0

    istring = qrm_str_tolower(string)

    if(index(istring,'qrm_amalgth') .eq. 1) then
       qrm_spmat%rcntl(qrm_amalgth_) = rval
    else if(index(istring,'qrm_rweight') .eq. 1) then
       qrm_spmat%rcntl(qrm_rweight_) = rval
    else if(index(istring,'qrm_mem_relax') .eq. 1) then
       qrm_spmat%rcntl(qrm_mem_relax_) = rval
    else
       err = 23
       call qrm_error_print(err, name, aed=string)
       goto 9999
    end if

9999 continue

    if(present(info)) info = err
    return

  end subroutine _qrm_psetr


  !> @brief Gets the values of an integer control parameter.
  !!        This is the dual of the @link ::_qrm_pseti @endlink
  !!        routine; the parameters and accepted values are the same.
  subroutine _qrm_pgeti(qrm_spmat, string, ival, info)
    use qrm_common_mod
    use qrm_string_mod
    use qrm_error_mod
    implicit none

    type(_qrm_spmat_type)           :: qrm_spmat
    character(len=*)                :: string
    integer                         :: ival
    integer, optional               :: info

    character(len=len(string))      :: istring
    integer(kind=8)                 :: iival

    ! error management
    integer                         :: err
    character(len=*), parameter     :: name='_qrm_pgeti'

    call _qrm_pgetii(qrm_spmat, string, iival, info)

    ival = iival

    return

  end subroutine _qrm_pgeti

  !> @brief Gets the values of an integer control parameter.
  !!        This is the dual of the @link ::_qrm_pseti @endlink
  !!        routine; the parameters and accepted values are the same.
  subroutine _qrm_pgetii(qrm_spmat, string, ival, info)
    use qrm_common_mod
    use qrm_string_mod
    use qrm_error_mod
    implicit none

    type(_qrm_spmat_type)           :: qrm_spmat
    character(len=*  )              :: string
    integer(kind=8)                 :: ival
    integer, optional               :: info

    character(len=len(string))      :: istring

    ! error management
    integer                         :: err
    character(len=*), parameter     :: name='_qrm_pgetii'

    err = 0

    istring = qrm_str_tolower(string)

    if(index(istring,'qrm_ordering') .eq. 1) then
       ival = qrm_spmat%icntl(qrm_ordering_)
    else if (index(istring,'qrm_minamalg') .eq. 1) then
       ival = qrm_spmat%icntl(qrm_minamalg_)
    else if (index(istring,'qrm_nb') .eq. 1) then
       ival = qrm_spmat%icntl(qrm_nb_)
    else if (index(istring,'qrm_mb') .eq. 1) then
       ival = qrm_spmat%icntl(qrm_mb_)
    else if (index(istring,'qrm_ib') .eq. 1) then
       ival = qrm_spmat%icntl(qrm_ib_)
    else if (index(istring,'qrm_bh') .eq. 1) then
       ival = qrm_spmat%icntl(qrm_bh_)
    else if (index(istring,'qrm_rhsnb') .eq. 1) then
       ival = qrm_spmat%icntl(qrm_rhsnb_)
    else if (index(istring,'qrm_keeph') .eq. 1) then
       ival = qrm_spmat%icntl(qrm_keeph_)
    else if (index(istring,'qrm_sing') .eq. 1) then
       ival = qrm_spmat%icntl(qrm_sing_)
    else if (index(istring,'qrm_e_nnz_r') .eq. 1) then
       ival = qrm_spmat%gstats(qrm_e_nnz_r_)
    else if (index(istring,'qrm_e_nnz_h') .eq. 1) then
       ival = qrm_spmat%gstats(qrm_e_nnz_h_)
    else if (index(istring,'qrm_e_facto_flops') .eq. 1) then
       ival = qrm_spmat%gstats(qrm_e_facto_flops_)
    else if (index(istring,'qrm_nnz_r') .eq. 1) then
       ival = qrm_spmat%gstats(qrm_nnz_r_)
    else if (index(istring,'qrm_nnz_h') .eq. 1) then
       ival = qrm_spmat%gstats(qrm_nnz_h_)
    else if (index(istring,'qrm_facto_flops') .eq. 1) then
       ival = qrm_spmat%gstats(qrm_facto_flops_)
    else
       err = 1
       call qrm_error_print(err, name, aed=string)
       goto 9999
    end if

9999 continue

    if(present(info)) info = err
    return

  end subroutine _qrm_pgetii



  !> @brief Gets the values of a real control parameter.
  !!        This is the dual of the @link ::_qrm_psetr @endlink
  !!        routine; the parameters and accepted values are the same.
  !!
  subroutine _qrm_pgetr(qrm_spmat, string, rval, info)

    use qrm_common_mod
    use qrm_string_mod
    use qrm_error_mod
    implicit none

    type(_qrm_spmat_type)           :: qrm_spmat
    character(len=*)                :: string
    real(kind(1.d0))                :: rval
    integer, optional               :: info
    
    character(len=len(string))      :: istring

    ! error management
    integer                         :: err
    character(len=*), parameter     :: name='_qrm_pgetr'

    err = 0
    
    istring = qrm_str_tolower(string)

    if(index(istring,'qrm_amalgth') .eq. 1) then
       rval = qrm_spmat%rcntl(qrm_amalgth_)
    else if(index(istring,'qrm_mem_relax') .eq. 1) then
       rval = qrm_spmat%rcntl(qrm_mem_relax_)
    else
       err = 1
       call qrm_error_print(err, name, aed=string)
       goto 9999
    end if

9999 continue

    if(present(info)) info = err
    return

  end subroutine _qrm_pgetr


  !> @brief Check the compatibility and correctness of icntl and rcntl
  !! parameters
  subroutine _qrm_check_spmat(qrm_spmat, op, info)
    use qrm_common_mod
    use qrm_string_mod
    use qrm_error_mod
    implicit none

    type(_qrm_spmat_type)           :: qrm_spmat
    integer, optional               :: op
    integer, optional               :: info

    integer                         :: iop

    ! error management
    integer                         :: err
    character(len=*), parameter     :: name='_qrm_check_spmat'

    err = 0

    if(present(op)) then
       iop = op
    else
       iop = qrm_allop_
    end if
    
    if((qrm_spmat%m .lt. 0) .or. (qrm_spmat%n .lt. 0) .or. &
         & (qrm_spmat%nz .lt. 0) .or. &
         & (qrm_spmat%nz .gt. (int(qrm_spmat%n,kind=8)*int(qrm_spmat%m,kind=8)))) then
       err = 29
       call qrm_error_print(err, name,ied=(/qrm_spmat%m,qrm_spmat%n,qrm_spmat%nz/))
       goto 9999
    end if


    if((iop.eq.qrm_allop_) .or. (iop.eq.qrm_analyse_)) then
       
       ! all the potential cases of conflict with the orderings
       select case(qrm_spmat%icntl(qrm_ordering_))
       case(:-1,6:)
          err = 9
          call qrm_error_print(err, name,ied=(/qrm_spmat%icntl(qrm_ordering_),0,0,0,0/))
          goto 9999
       case (qrm_given_)
          if(qrm_spmat%icntl(qrm_sing_) .eq. qrm_yes_) then
             err = 27
             call qrm_error_print(err, name,ied=(/qrm_ordering_,qrm_sing_/))
             goto 9999
          end if
       end select
       
       ! all the potential cases of conflict with the orderings
       select case(qrm_spmat%icntl(qrm_nb_))
       case(:-1)
          err = 28
          call qrm_error_print(err, name, ied=(/qrm_spmat%icntl(qrm_mb_),&
               & qrm_spmat%icntl(qrm_nb_),qrm_spmat%icntl(qrm_ib_)/))
          goto 9999
       case default
          if(qrm_spmat%icntl(qrm_nb_) .lt. qrm_spmat%icntl(qrm_ib_)) then
             err = 27
             call qrm_error_print(err, name, ied=(/qrm_nb_,qrm_ib_/))
             goto 9999
          end if
       end select

       if(    (mod(qrm_spmat%icntl(qrm_mb_),qrm_spmat%icntl(qrm_nb_)).gt.0) .or. &
            & (mod(qrm_spmat%icntl(qrm_nb_),qrm_spmat%icntl(qrm_ib_)).gt.0) .or. &
            & ( (qrm_spmat%icntl(qrm_mb_).gt.0) .and. (qrm_spmat%icntl(qrm_mb_) .lt. qrm_spmat%icntl(qrm_nb_)))) then
          err = 28
          call qrm_error_print(err, name, ied=(/qrm_spmat%icntl(qrm_mb_),&
               & qrm_spmat%icntl(qrm_nb_),qrm_spmat%icntl(qrm_ib_)/))
          goto 9999
       end if

       select case(qrm_spmat%icntl(qrm_ib_))
       case(:-1)
          err = 28
          call qrm_error_print(err, name, ied=(/qrm_spmat%icntl(qrm_ib_)/))
          goto 9999
       end select
    end if

9999 continue

    if(present(info)) info = err
    return

  end subroutine _qrm_check_spmat
  

  subroutine _qrm_sync_spmat(qrm_mat)
    use qrm_common_mod
    use qrm_error_mod
    use _qrm_fdata_mod
    use qrm_mem_mod
#if defined(have_starpu)
    use starpu_f_mod
#endif
    implicit none

    type(_qrm_spmat_type), target  :: qrm_mat

    type(_qrm_front_type), pointer :: front
    type(qrm_adata_type), pointer  :: adata
    type(_qrm_fdata_type), pointer :: fdata

    integer :: inode, node

    adata => qrm_mat%adata
    fdata => qrm_mat%fdata
    
#if defined(have_starpu)
    if(allocated(qrm_mat%adata)) then
       if(c_associated(adata%hdl)) then
          call starpu_f_data_acquire_read(adata%hdl)
          call starpu_f_data_release(adata%hdl)
       end if
    end if

    if(allocated(qrm_mat%fdata)) then
       if(allocated(fdata%front_list)) then
          do inode = 1, adata%nnodes
             node = adata%torder(inode)
             if(adata%small(node) .ge. 0) then
                front => fdata%front_list(node)
                if(c_associated(front%sym_handle)) then
                   call starpu_f_data_acquire_read(front%sym_handle)
                   call starpu_f_data_release(front%sym_handle)
                end if
             end if
          end do
       end if
    end if
#endif
    return

  end subroutine _qrm_sync_spmat
  
end module _qrm_spmat_mod
