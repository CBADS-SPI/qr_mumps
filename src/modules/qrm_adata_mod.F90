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
!> @file qrm_adata_mod.F90
!! This file contains the module that holds the main data type for the analysis phase.
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


#include "qrm_common.h"


!> @brief This module contains the definition of the analysis data type.



module qrm_adata_mod

#if defined (have_starpu)
  use starpu_f_mod
#endif
  
  !> @brief The main data type for the analysis phase.

  !> This data type is meant to store all the results of the anaysis
  !! phase (permutations, assembly tree etc.). Moreover it is
  !! not typed which means that it will be exactly the same for dp/sp
  !! real or complex matrices. However, because of its untyped nature,
  !! this data type cannot be included into qrm_analysis_mod which, instead
  !! is typed due to the fact that many routines there take as input
  !! argument a qrm_spmat_type.

  type qrm_adata_type
     !> An integer array of size n (column dimension of the problem matrix)
     !! containing the column permutation resulting from fill-in
     !! minimization orderings (COLAMD etc). It has to be a pointer to
     !! accomodate the case where the user wants to provide his own
     !! permutation.
     integer, allocatable            :: cperm(:)
     !> This is an integer array containing a row permutation computed in order
     !! to have the original matrix in a "stair" format. See @link _qrm_rowperm_ @endlink
     integer, allocatable            :: rperm(:)
     !> This integer array holds pointers to the beginning of nodes in the
     !! column permutation. The variables contained in node i are
     !! cperm( cp_ptr(i) : cp_ptr(i+1)-1 ). It is of size nnodes+1
     integer, allocatable            :: cp_ptr(:)
     !> rc(i) contains the rowcount for the rows in node number i. This array
     !! is of size nnodes
     integer, allocatable            :: rc(:)
     !> An integer array of size nnodes containing the assembly tree.
     !! parent(i)=k means that node k is the father of node i
     integer, allocatable            :: parent(:)
     !> For each node, it contains the list of its children
     integer, allocatable            :: child(:)
     !> Pointers to the list of children. child(childptr(i),...,childptr(i+1)-1)
     !! contains all the children of node i
     integer, allocatable            :: childptr(:)
     !> nfrows(i)=k means that the frontal matrix related to node i has k rows
     !! (the number of columns is rc(i))
     integer, allocatable            :: nfrows(:)
     !> defines the staircase structure of the original matrix. Assuming A(rperm,cperm)
     !! @verbatim
     !!            |x       |
     !!            |x       |     in this case stair would be:
     !!            | x      |     stair=(/ 2, 5, 7, 7, 9 /)
     !!            | x      |
     !!            | x      |
     !!            |  xx    |
     !!            |  xx    |
     !!            |    x   |
     !!            |    x   |
     !! @endverbatim         
     integer, allocatable            :: stair(:)
     !> An array of size nnodes which flags (with 1 whereas everything is =0) the roots
     !! of subtrees that are treated sequentially in the numerical factorization
     integer, allocatable            :: small(:) 
     !> Contains the list of columns indices for each front
     integer, allocatable            :: fcol(:) 
     !> Contains pointers to the list of column indices fcol. Specifically, the list of
     !! column indices for front i is fcol(fcol_ptr(i):fcol_ptr(i+1)-1)
     integer, allocatable            :: fcol_ptr(:)
     !> A list of leaf nodes where to start the numerical factorization from. A leaf is
     !! defined as a node which only has small children
     integer, allocatable            :: leaves(:) 
     !> The tree traversal order
     integer, allocatable            :: torder(:) 
     !> arrays containing the memory needed upon activation and freed upon clean
     integer(kind=8), allocatable    :: asize(:), csize(:)
     !> array containing the number of flops per node
     integer(kind=8), allocatable    :: weight(:)
     !> The number of leaves present in the tree
     integer                         :: nleaves=0
     !> The number of nodes in the elimination tree, i.e., the number of frontal
     !! matrices.
     integer                         :: nnodes=0
     !> it is set to .true. if the analysis is done (with success), to .false. otherwise
     logical                         :: ok=.false.

#if defined (have_starpu)
     type(c_ptr)                     :: hdl
#endif
  end type qrm_adata_type

contains

  !> @brief initializes a @link qrm_adata_type @endlink instance
  !! @param[inout]  adata_in  the instance to be initialized
  subroutine qrm_adata_init(qrm_adata, info, register)

    use qrm_error_mod
    implicit none

    type(qrm_adata_type), allocatable :: qrm_adata
    integer, optional                 :: info
    logical, optional                 :: register
    
    logical                           :: iregister

    ! error management
    integer                           :: err
    character(len=*), parameter       :: name='qrm_adata_init'

    err = 0

    allocate(qrm_adata, stat=err)
#if defined (have_starpu)
    if(err.eq.0) then
       if(present(register)) then
          iregister = register
       else
          iregister = .true.
       end if

       if(iregister) then
          call starpu_f_void_data_register(qrm_adata%hdl)
       else
          qrm_adata%hdl = c_null_ptr
       end if
    end if
#endif    

    if(present(info)) info = err
    return
    
  end subroutine qrm_adata_init
  
  
  !> @brief make s a copy of an @link qrm_adata_type @endlink instance
  !! @param[in]  adata_in  the instance to be copied
  !! @param[out] adata_out the produced copy
  subroutine qrm_adata_copy(qrm_adata_in, qrm_adata_out, info)

    use qrm_mem_mod
    use qrm_error_mod
    implicit none

    type(qrm_adata_type), intent(in)  :: qrm_adata_in
    type(qrm_adata_type), intent(out) :: qrm_adata_out
    integer, optional                 :: info
    
    ! error management
    integer                           :: err
    character(len=*), parameter       :: name='qrm_adata_copy'

    err = 0

    if(err.eq.0) call qrm_alloc(qrm_adata_out%cperm,    qrm_size(qrm_adata_in%cperm),    err)
    if(err.eq.0) call qrm_alloc(qrm_adata_out%rperm,    qrm_size(qrm_adata_in%rperm),    err)
    if(err.eq.0) call qrm_alloc(qrm_adata_out%cp_ptr,   qrm_size(qrm_adata_in%cp_ptr),   err)
    if(err.eq.0) call qrm_alloc(qrm_adata_out%rc,       qrm_size(qrm_adata_in%rc),       err)
    if(err.eq.0) call qrm_alloc(qrm_adata_out%parent,   qrm_size(qrm_adata_in%parent),   err)
    if(err.eq.0) call qrm_alloc(qrm_adata_out%child,    qrm_size(qrm_adata_in%child),    err)
    if(err.eq.0) call qrm_alloc(qrm_adata_out%childptr, qrm_size(qrm_adata_in%childptr), err)
    if(err.eq.0) call qrm_alloc(qrm_adata_out%nfrows,   qrm_size(qrm_adata_in%nfrows),   err)
    if(err.eq.0) call qrm_alloc(qrm_adata_out%stair,    qrm_size(qrm_adata_in%stair),    err)
    if(err.eq.0) call qrm_alloc(qrm_adata_out%leaves,   qrm_size(qrm_adata_in%leaves),   err)
    if(err.eq.0) call qrm_alloc(qrm_adata_out%fcol,     qrm_size(qrm_adata_in%fcol),     err)
    if(err.eq.0) call qrm_alloc(qrm_adata_out%fcol_ptr, qrm_size(qrm_adata_in%fcol_ptr), err)
    if(err.eq.0) call qrm_alloc(qrm_adata_out%small,    qrm_size(qrm_adata_in%small),    err)
    if(err.eq.0) call qrm_alloc(qrm_adata_out%asize,    qrm_size(qrm_adata_in%asize),    err)
    if(err.eq.0) call qrm_alloc(qrm_adata_out%csize,    qrm_size(qrm_adata_in%csize),    err)
    __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)    

    if(allocated(qrm_adata_in%cperm))qrm_adata_out%cperm       = qrm_adata_in%cperm 
    if(allocated(qrm_adata_in%rperm))qrm_adata_out%rperm       = qrm_adata_in%rperm
    if(allocated(qrm_adata_in%cp_ptr))qrm_adata_out%cp_ptr     = qrm_adata_in%cp_ptr
    if(allocated(qrm_adata_in%rc))qrm_adata_out%rc             = qrm_adata_in%rc
    if(allocated(qrm_adata_in%parent))qrm_adata_out%parent     = qrm_adata_in%parent
    if(allocated(qrm_adata_in%child))qrm_adata_out%child       = qrm_adata_in%child
    if(allocated(qrm_adata_in%childptr))qrm_adata_out%childptr = qrm_adata_in%childptr
    if(allocated(qrm_adata_in%nfrows))qrm_adata_out%nfrows     = qrm_adata_in%nfrows
    if(allocated(qrm_adata_in%stair))qrm_adata_out%stair       = qrm_adata_in%stair
    if(allocated(qrm_adata_in%leaves))qrm_adata_out%leaves     = qrm_adata_in%leaves
    if(allocated(qrm_adata_in%fcol))qrm_adata_out%fcol         = qrm_adata_in%fcol
    if(allocated(qrm_adata_in%fcol_ptr))qrm_adata_out%fcol_ptr = qrm_adata_in%fcol_ptr
    if(allocated(qrm_adata_in%small))qrm_adata_out%small       = qrm_adata_in%small
    if(allocated(qrm_adata_in%asize))qrm_adata_out%asize       = qrm_adata_in%asize
    if(allocated(qrm_adata_in%csize))qrm_adata_out%csize       = qrm_adata_in%csize

    qrm_adata_out%nnodes     = qrm_adata_in%nnodes
    qrm_adata_out%nleaves    = qrm_adata_in%nleaves
    qrm_adata_out%ok         = qrm_adata_in%ok

9999 continue 
    if(present(info)) info = err
    return

  end subroutine qrm_adata_copy


  !> @brief make s a move of an @link qrm_adata_type @endlink instance
  !! @param[in]  adata_in  the instance to be moved
  !! @param[out] adata_out the produced copy
  subroutine qrm_adata_move(qrm_adata_in, qrm_adata_out, info)

    use qrm_mem_mod
    use qrm_error_mod
    implicit none

    type(qrm_adata_type), intent(inout)  :: qrm_adata_in
    type(qrm_adata_type), intent(inout)  :: qrm_adata_out
    integer, optional                    :: info

    ! error management
    integer                              :: err
    character(len=*), parameter          :: name='qrm_adata_move'

    err = 0

    call move_alloc(qrm_adata_in%cperm   ,  qrm_adata_out%cperm     ) 
    call move_alloc(qrm_adata_in%rperm   ,  qrm_adata_out%rperm     ) 
    call move_alloc(qrm_adata_in%cp_ptr  ,  qrm_adata_out%cp_ptr    ) 
    call move_alloc(qrm_adata_in%rc      ,  qrm_adata_out%rc        ) 
    call move_alloc(qrm_adata_in%parent  ,  qrm_adata_out%parent    ) 
    call move_alloc(qrm_adata_in%child   ,  qrm_adata_out%child     ) 
    call move_alloc(qrm_adata_in%childptr,  qrm_adata_out%childptr  ) 
    call move_alloc(qrm_adata_in%nfrows  ,  qrm_adata_out%nfrows    ) 
    call move_alloc(qrm_adata_in%stair   ,  qrm_adata_out%stair     ) 
    call move_alloc(qrm_adata_in%leaves  ,  qrm_adata_out%leaves    ) 
    call move_alloc(qrm_adata_in%fcol    ,  qrm_adata_out%fcol      ) 
    call move_alloc(qrm_adata_in%fcol_ptr,  qrm_adata_out%fcol_ptr  ) 
    call move_alloc(qrm_adata_in%small   ,  qrm_adata_out%small     ) 
    call move_alloc(qrm_adata_in%torder  ,  qrm_adata_out%torder    ) 
    call move_alloc(qrm_adata_in%asize   ,  qrm_adata_out%asize     ) 
    call move_alloc(qrm_adata_in%csize   ,  qrm_adata_out%csize     ) 
    call move_alloc(qrm_adata_in%weight  ,  qrm_adata_out%weight    ) 

    qrm_adata_out%nnodes     = qrm_adata_in%nnodes
    qrm_adata_out%nleaves    = qrm_adata_in%nleaves
    qrm_adata_out%ok         = qrm_adata_in%ok

    if(present(info)) info = err
    return

  end subroutine qrm_adata_move


  !> @brief Frees an @link qrm_adata_type @endlink instance
  !! @param[in,out] adata the instance to be freed
  subroutine qrm_adata_destroy(qrm_adata, info)

    use qrm_mem_mod
    use qrm_error_mod
    implicit none

    type(qrm_adata_type), allocatable :: qrm_adata
    integer, optional                 :: info

    ! error management
    integer                           :: err
    character(len=*), parameter       :: name='qrm_adata_destroy'

    err = 0

    if(.not.allocated(qrm_adata)) goto 9999
    
    call qrm_adata_cleanup(qrm_adata, err)
    __QRM_INFO_CHECK(err, name,'qrm_cleanup',9999)

#if defined (have_starpu)
    if(c_associated(qrm_adata%hdl)) then
       call starpu_f_data_acquire_read(qrm_adata%hdl)
       call starpu_f_data_release(qrm_adata%hdl)
       call starpu_f_data_unregister(qrm_adata%hdl)
    end if
#endif
    deallocate(qrm_adata)
    
9999 continue 
    if(present(info)) info = err
    return

  end subroutine qrm_adata_destroy


  !> @brief cleans up an @link qrm_adata_type @endlink instance
  !! @param[in,out] adata the instance to be freed
  subroutine qrm_adata_cleanup(qrm_adata, info)
    use qrm_mem_mod
    use qrm_error_mod
    implicit none

    type(qrm_adata_type)            :: qrm_adata
    integer, optional               :: info
  
    ! error management
    integer                         :: err
    character(len=*), parameter     :: name='qrm_adata_cleanup'

    err = 0

    if(err.eq.0) call qrm_dealloc(qrm_adata%cperm,    err)
    if(err.eq.0) call qrm_dealloc(qrm_adata%rperm,    err)
    if(err.eq.0) call qrm_dealloc(qrm_adata%cp_ptr,   err)
    if(err.eq.0) call qrm_dealloc(qrm_adata%rc,       err)
    if(err.eq.0) call qrm_dealloc(qrm_adata%parent,   err)
    if(err.eq.0) call qrm_dealloc(qrm_adata%nfrows,   err)
    if(err.eq.0) call qrm_dealloc(qrm_adata%stair,    err)
    if(err.eq.0) call qrm_dealloc(qrm_adata%small,    err)
    if(err.eq.0) call qrm_dealloc(qrm_adata%childptr, err)
    if(err.eq.0) call qrm_dealloc(qrm_adata%child,    err)
    if(err.eq.0) call qrm_dealloc(qrm_adata%leaves,   err)
    if(err.eq.0) call qrm_dealloc(qrm_adata%fcol,     err)
    if(err.eq.0) call qrm_dealloc(qrm_adata%fcol_ptr, err)
    if(err.eq.0) call qrm_dealloc(qrm_adata%asize,    err)
    if(err.eq.0) call qrm_dealloc(qrm_adata%csize,    err)
    if(err.eq.0) call qrm_dealloc(qrm_adata%torder,   err)
    if(err.eq.0) call qrm_dealloc(qrm_adata%weight,   err)
    __QRM_INFO_CHECK(err, name,'qrm_dealloc',9999)

    qrm_adata%nleaves = 0
    qrm_adata%nnodes  = 0
    qrm_adata%ok = .false.

9999 continue 
    if(present(info)) info = err
    return

  end subroutine qrm_adata_cleanup
  
    
end module qrm_adata_mod
