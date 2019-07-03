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
!> @file qrm_fdata_mod.F90
!! This file contains the module which holds the data types used during the factorization.
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief This module contains the definition of all the data related to the
!! factorization phase. 
!! 
module _qrm_fdata_mod
#if defined (have_starpu)
  use starpu_f_mod
#endif
  use qrm_pthread_mod
  use qrm_memhandling_mod
  
  !> @brief This type defines a block-column (a simple 2D array
  !> containing the coefficients of a block-column in the front)
  type _qrm_bc_type
     !> the coefficients
     _qrm_data, allocatable          :: c(:,:)
#if defined (have_starpu)
     type(c_ptr)    :: hdl = c_null_ptr
#endif
  end type _qrm_bc_type
  
  !> @brief This type defines a data structure containing all the data related to a front.
  type _qrm_front_type
     ! amount of flop performed on the front 
     real(kind(1.d0))                :: flops = 0.d0
     ! number of task   
     integer(kind=8)                 :: ntask = 0

     !> the handle representing symbolic data in the front
#if defined (have_starpu)
     type(c_ptr)                     :: sym_handle
#endif
     !> the front number
     integer                         :: num=0
     !> the number of rows in the front
     integer                         :: m=0
     !> the number of columns in the front
     integer                         :: n=0
     !> the number of pivots (fully assembled variables in the front
     integer                         :: npiv=0
     !> the list of row indices in the front
     integer, allocatable            :: rows(:)
     !> the list of column indices in the front
     integer, allocatable            :: cols(:)
     !> The pointer to the beginning of rows of coefficients from
     !> the original matrix. 
     integer, allocatable            :: aiptr(:)
     !> The columns indices of coefficients from
     !> the original matrix. 
     integer, allocatable            :: ajcn(:)
     !> The values of coefficients from original matrix. These
     !! coefficients will be assembled into
     !! the front at the moment of its activation
     _qrm_data, allocatable          :: aval(:)
     !> The initial data from the problem matrix is stored in CSR format. This
     !  integer contains the number of rows in this front from the original matrix
     integer                         :: anrows
     !> The row mapping for the rows from
     !> the original matrix. 
     integer, allocatable            :: arowmap(:)
     !> this integer array of size front%n contains a mapping between the front's
     !! column indices and father front's column indices, i.e. colmap(i)=k means
     !! that column i of front goes into k-th column of its parent
     integer, allocatable            :: colmap(:)
     !> this integer array is of size equal to the number of rows in
     !! the contribution block produced by this front and contains a
     !! mapping between the rows of the CB and the front of the
     !! parent, i.e., rowmap(i)=k means that i-th row of this
     !! contribution block goes into row k of the front of the father.
     integer, allocatable            :: rowmap(:)
     !> this array of size front%n will ultimately define the staircase
     !! structure of the front. stair(i)=k means that the lowest
     !! element in column i of the front is in row k
     integer, allocatable            :: stair(:)
     !> an array containing all the blocks/tiles
     type(_qrm_bc_type), allocatable :: bc(:,:)
     !> after the front is factorized, the corresponding elemnts of R will be
     !! copied into this array
     type(_qrm_bc_type), allocatable :: r(:,:)
     !> after the front is factorized, the corresponding elemnts of H will be
     !! copied into this array
     type(_qrm_bc_type), allocatable :: h(:,:)
     !> these blocks hold the T matrices resulting from each block factorization. 
     type(_qrm_bc_type), allocatable :: t(:,:), t2(:,:)
     !> the block sizes for this front
     integer                         :: mb, nb, ib, bh
     !> the number of column and row-blocks in this front
     integer                         :: nc, nr
     !> the number of block-panels to be factorized
     integer                         :: np
     !> the number of eliminations to be performed on the front.
     !!  This basically corresponds to min(front%m,front%n)
     integer                         :: ne 
     !> the number of entries in R and H in this front
     integer(kind=8)                 :: rsize, hsize

  end type _qrm_front_type



  !> @brief The data structure meant to store all the results of the
  !! factorization phase.
  type _qrm_fdata_type
     !> an integer containing the number of fronts assigned to the process
     integer                             :: nfronts=0
     ! !> just a counter to count the number of fronts done during the
     ! !> factorization. It is used to stop the factorization operation
     ! integer                             :: done=0
     !> an array of type _qrm_front_type containing the list of fronts
     !! assigned to the process
     type(_qrm_front_type), allocatable  :: front_list(:)
     !> it is set to .true. if the factorization is done (with success), to .false. otherwise
     logical                             :: ok=.false.
     !> workspace various factorization tasks
     type(_qrm_bc_type)                  :: work
     !> data structure to keep track of mem consumption in memaware execution
     type(qrm_ma_type)                   :: ma
  end type _qrm_fdata_type

  interface qrm_fdata_destroy
     module procedure _qrm_fdata_destroy
  end interface qrm_fdata_destroy

  interface qrm_fdata_cleanup
     module procedure _qrm_fdata_cleanup
  end interface qrm_fdata_cleanup

  interface qrm_fdata_init
     module procedure _qrm_fdata_init
  end interface qrm_fdata_init



contains

  !> @brief Frees a qrm_front_type instance
  !!
  !! @param[in,out] qrm_front the data to be freed
  !!
  subroutine _qrm_front_destroy(qrm_front, info)
    use qrm_mem_mod
    use qrm_error_mod
#if defined (have_starpu)
    use starpu_f_mod
#endif
    implicit none

    type(_qrm_front_type)       :: qrm_front
    integer, optional           :: info
    
    integer                     :: i, j
    ! error management
    integer                     :: err
    character(len=*), parameter :: name='qrm_front_destroy'
    
    err = 0
    
    if(err.eq.0) call qrm_dealloc(qrm_front%aiptr,   err)
    if(err.eq.0) call qrm_dealloc(qrm_front%ajcn,    err)
    if(err.eq.0) call qrm_dealloc(qrm_front%aval,    err)
    if(err.eq.0) call qrm_dealloc(qrm_front%rows,    err)
    if(err.eq.0) call qrm_dealloc(qrm_front%cols,    err)
    if(err.eq.0) call qrm_dealloc(qrm_front%arowmap, err)
    if(err.eq.0) call qrm_dealloc(qrm_front%rowmap,  err)
    if(err.eq.0) call qrm_dealloc(qrm_front%colmap,  err)
    if(err.eq.0) call qrm_dealloc(qrm_front%stair,   err)
    __QRM_INFO_CHECK(err, name, 'qrm_dealloc', 9999)
    
    if(allocated(qrm_front%bc)) then
       do i=1, size(qrm_front%bc,1)
          do j=1, size(qrm_front%bc,2)
             call qrm_dealloc(qrm_front%bc(i,j)%c, err)
             __QRM_INFO_CHECK(err, name, 'qrm_dealloc', 9999)
          end do
       end do
       deallocate(qrm_front%bc)
    end if
    if(allocated(qrm_front%t)) then
       do i=1, size(qrm_front%t,1)
          do j=1, size(qrm_front%t,2)
             call qrm_dealloc(qrm_front%t(i,j)%c, err)
             __QRM_INFO_CHECK(err, name, 'qrm_dealloc', 9999)
          end do
       end do
       deallocate(qrm_front%t)
    end if
    if(allocated(qrm_front%t2)) then
       do i=1, size(qrm_front%t2,1)
          do j=1, size(qrm_front%t2,2)
             call qrm_dealloc(qrm_front%t2(i,j)%c, err)
             __QRM_INFO_CHECK(err, name, 'qrm_dealloc', 9999)
          end do
       end do
       deallocate(qrm_front%t2)
    end if
    if(allocated(qrm_front%r)) then
       do i=1, size(qrm_front%r,1)
          do j=1, size(qrm_front%r,2)
             call qrm_dealloc(qrm_front%r(i,j)%c, err)
             __QRM_INFO_CHECK(err, name, 'qrm_dealloc', 9999)
          end do
       end do
       deallocate(qrm_front%r)
    end if
    if(allocated(qrm_front%h)) then
       do i=1, size(qrm_front%h,1)
          do j=1, size(qrm_front%h,2)
             call qrm_dealloc(qrm_front%h(i,j)%c, err)
             __QRM_INFO_CHECK(err, name, 'qrm_dealloc', 9999)
          end do
       end do
       deallocate(qrm_front%h)
    end if

#if defined (have_starpu)
    if(c_associated(qrm_front%sym_handle)) call starpu_f_data_unregister(qrm_front%sym_handle)
#endif
    
    qrm_front%m = 0
    qrm_front%n = 0
    
    
9999 continue
    if(present(info)) info = err
    return
    
  end subroutine _qrm_front_destroy
  


  !> @brief Inits a @link _qrm_fdata_type @endlink instance
  !! 
  !! @param[in,out] qrm_fdata the instace to be inited
  !! 
  subroutine _qrm_fdata_init(qrm_fdata, info)
    use qrm_mem_mod
    use qrm_error_mod
    implicit none
    
    type(_qrm_fdata_type), allocatable :: qrm_fdata
    integer, optional                  :: info
                                        
    integer                            :: i
                                        
    ! error management                  
    integer                            :: err
    character(len=*), parameter        :: name='qrm_fdata_init'
    
    err = 0

    allocate(qrm_fdata, stat=err)
    if(err.eq.0) then
       qrm_fdata%nfronts=0
       qrm_fdata%ok = .false.
    end if
    
    if(present(info)) info = err
    return
    
  end subroutine _qrm_fdata_init
  
    
  !> @brief Destroys a @link _qrm_fdata_type @endlink instance
  !! 
  !! @param[in,out] qrm_fdata the instace to be freed
  !! 
  subroutine _qrm_fdata_destroy(qrm_fdata, info)
    use qrm_mem_mod
    use qrm_error_mod
    implicit none
    
    type(_qrm_fdata_type), allocatable :: qrm_fdata
    integer, optional                  :: info
                                        
    integer                            :: i
                                        
    ! error management                  
    integer                            :: err
    character(len=*), parameter        :: name='qrm_fdata_destroy'
    
    err = 0
    
    if(allocated(qrm_fdata)) then
       call _qrm_fdata_cleanup(qrm_fdata, err)
       deallocate(qrm_fdata)
    end if
    
9999 continue
    if(present(info)) info = err
    return
    
  end subroutine _qrm_fdata_destroy
  
  !> @brief Destroys a @link _qrm_fdata_type @endlink instance
  !! 
  !! @param[in,out] qrm_fdata the instace to be freed
  !! 
  subroutine _qrm_fdata_cleanup(qrm_fdata, info)
    use qrm_mem_mod
    use qrm_error_mod
    implicit none
    
    type(_qrm_fdata_type)           :: qrm_fdata
    integer, optional               :: info

    integer                         :: i
    
    ! error management
    integer                         :: err
    character(len=*), parameter     :: name='qrm_fdata_cleanup'
    
    err = 0
    
    if(allocated(qrm_fdata%front_list)) then
       do i=1, qrm_fdata%nfronts
          call _qrm_front_destroy(qrm_fdata%front_list(i), err)
       end do
       deallocate(qrm_fdata%front_list)
    end if
    __QRM_INFO_CHECK(err, name,'qrm_front_destroy',9999)

    call qrm_facto_mem_finalize(qrm_fdata%ma)
    
    ! qrm_fdata%done    = 0
    qrm_fdata%nfronts = 0
    qrm_fdata%ok      = .false.

9999 continue
    if(present(info)) info = err
    return
    
  end subroutine _qrm_fdata_cleanup
  


end module _qrm_fdata_mod
