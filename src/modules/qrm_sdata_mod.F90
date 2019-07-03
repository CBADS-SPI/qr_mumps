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
!> @file qrm_sdata_mod.F90
!! This file contains a module that defines the data type for storing
!! the results of the solve phase (not used at the moment)
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


!> @brief This module contains the definition of all the data related to the
!! solution phase. 
!! 
module _qrm_sdata_mod
  use _qrm_fdata_mod     ! TODO add deps in the makefile
#if defined (have_starpu)
  use starpu_f_mod
#endif


  ! Type: _qrm_rhs_type
  !
  ! This type defines a data structure containing a rhs block
  !
  type _qrm_rhs_type
     _qrm_data, pointer              :: p(:,:)
#if defined (have_starpu)
     type(c_ptr)                     :: handle
#endif
     type(_qrm_bc_type), allocatable :: front_rhs(:)
     type(_qrm_bc_type)              :: work
     type(_qrm_fdata_type), pointer  :: fdata_pointer
  end type _qrm_rhs_type

  interface qrm_rhs_init
     module procedure _qrm_rhs_init2d, _qrm_rhs_init1d
  end interface qrm_rhs_init

  interface _qrm_rhs_init
     module procedure _qrm_rhs_init2d, _qrm_rhs_init1d
  end interface _qrm_rhs_init

  interface qrm_rhs_destroy
     module procedure _qrm_rhs_destroy
  end interface qrm_rhs_destroy

  interface qrm_sync_rhs
     module procedure  _qrm_sync_rhs
  end interface qrm_sync_rhs

  interface qrm_sync_data
     module procedure :: _qrm_sync_rhs
  end interface qrm_sync_data  
  
contains
  
  subroutine _qrm_rhs_init2d(qrm_rhs, x)
#if defined (have_starpu)
    use starpu_f_mod
#endif
    use qrm_mem_mod
    implicit none

    type(_qrm_rhs_type) :: qrm_rhs
    _qrm_data, target   :: x(:,:)
    
    qrm_rhs%p => x

#if defined (have_starpu)
    call starpu_f_matrix_data_register(qrm_rhs%handle, 0, &
         & qrm_rhs%p, &
         & size(qrm_rhs%p,1), &
         & size(qrm_rhs%p,1), &
         & size(qrm_rhs%p,2))
#endif

    if(allocated(qrm_rhs%front_rhs)) deallocate(qrm_rhs%front_rhs)
    call qrm_dealloc(qrm_rhs%work%c)

    return
  end subroutine _qrm_rhs_init2d

  subroutine _qrm_rhs_init1d(qrm_rhs, x)
    use qrm_mem_mod
    use _qrm_utils_mod
    implicit none

    type(_qrm_rhs_type) :: qrm_rhs
    _qrm_data, target   :: x(:)
    _qrm_data, pointer  :: pnt(:,:)
    
    ! Turn the 1D pointer into a 2D one and call the 2D method below
    call _qrm_remap_pnt(x, pnt, size(x,1))    

    call _qrm_rhs_init2d(qrm_rhs, pnt)

  end subroutine _qrm_rhs_init1d

  
  
  
  subroutine _qrm_rhs_destroy(qrm_rhs)
#if defined (have_starpu)
    use starpu_f_mod
#endif
    use qrm_mem_mod
    use iso_c_binding
    implicit none

    type(_qrm_rhs_type) :: qrm_rhs

    integer             :: i

#if defined (have_starpu)
    call starpu_f_data_unregister(qrm_rhs%handle)
    qrm_rhs%handle = c_null_ptr
#endif

    if(allocated(qrm_rhs%front_rhs)) then
       do i=1, size(qrm_rhs%front_rhs,1)
#if defined (have_starpu)
          if(c_associated(qrm_rhs%front_rhs(i)%hdl)) then
             call starpu_f_data_unregister(qrm_rhs%front_rhs(i)%hdl)
             qrm_rhs%front_rhs(i)%hdl = c_null_ptr
          end if
#endif
          call qrm_dealloc(qrm_rhs%front_rhs(i)%c)
       end do
       deallocate(qrm_rhs%front_rhs)
    end if

    nullify(qrm_rhs%p)

    return
  end subroutine _qrm_rhs_destroy

  subroutine _qrm_sync_rhs(qrm_rhs)
    use qrm_common_mod
    use qrm_error_mod
#if defined(have_starpu)
    use starpu_f_mod
#endif
    implicit none
    
    type(_qrm_rhs_type), target  :: qrm_rhs
    
    integer :: inode, node
    
#if defined(have_starpu)
    if(allocated(qrm_rhs%front_rhs)) then
       do inode = 1, size(qrm_rhs%front_rhs,1)
          if(c_associated(qrm_rhs%front_rhs(inode)%hdl)) then
             call starpu_f_data_acquire_read(qrm_rhs%front_rhs(inode)%hdl)
             call starpu_f_data_release(qrm_rhs%front_rhs(inode)%hdl)
          end if
       end do
    end if
#endif
    return
    
  end subroutine _qrm_sync_rhs
  

end module _qrm_sdata_mod
