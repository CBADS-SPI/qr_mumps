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
!> @file qrm_do_colamd.F90
!! This file contains the routine that computes a COLAMD permutation of the input matrix
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


#include "qrm_common.h"


!> @brief This subroutine computes the fill reducing ordering using COLAMD

!> Please refer to:
!!
!! <em>"A column approximate minimum degree ordering algorithm"</em>,
!! T. A. Davis, J. R. Gilbert, S. Larimore, E. Ng, ACM Transactions on
!! Mathematical Software, vol 30, no. 3, Sept. 2004, pp. 353-376.
!! 
!! for the details of the reordering method.
!!
!! @param[in] graph  the graph associated to the matrix to be ordered.
!!
!! @param[out] cperm an integer array containing the new column order
!!
subroutine _qrm_do_colamd(graph, cperm, info)

  use _qrm_spmat_mod
  use qrm_error_mod
  use _qrm_analysis_mod, savesym => _qrm_do_colamd
  use qrm_mem_mod
  use iso_c_binding

  implicit none

  type(_qrm_spmat_type)       :: graph
  integer, target             :: cperm(:)
  integer, optional           :: info

  interface
     subroutine qrm_colamd(n_row, n_col, Alen, A, p, err) bind(c, name='qrm_colamd')
       use iso_c_binding
       integer(c_int), value :: n_row, n_col, Alen
       integer(c_int)        :: A(*), p(*), err
     end subroutine qrm_colamd
  end interface

  interface
     subroutine qrm_colamd_recommended(Alen, nnz, n_row, n_col) bind(c, name='qrm_colamd_recommended')
       use iso_c_binding
       integer(c_int), value :: nnz, n_row, n_col
       integer(c_int)        :: Alen
     end subroutine qrm_colamd_recommended
  end interface

  type(_qrm_spmat_type)           :: gcopy
  integer                         :: i, idx, cnt, tmp, alen
  ! error management
  integer                         :: err, err2
  character(len=*), parameter     :: name='qrm_do_colamd'

  err = 0; err2 = 0

  ! at this point we have to make a copy of the graph.
  ! this is a huge waste of mem but we don't have a choice
  ! since ccolamd destroys the graph which, instead, we want to
  ! save for successive computations
  
  ! compute the memory required by ccolamd (a lot) and allocate
  call qrm_colamd_recommended(alen, graph%nz, graph%m, graph%n)
  call qrm_alloc(gcopy%irn, alen, err)
  __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)
  call qrm_alloc(gcopy%jptr, graph%n+1)

  ! make the copy
  call _qrm_spmat_copy(graph, gcopy, values=.false., info=err)
  __QRM_INFO_CHECK(err, name,'qrm_spmat_copy',9999)

  ! ccolamd wants 0 based indices (argh!)
  gcopy%irn(1:gcopy%nz) = gcopy%irn(1:gcopy%nz)-1
  gcopy%jptr(1:gcopy%n+1) = gcopy%jptr(1:gcopy%n+1)-1

  ! call ccolamd
  call qrm_colamd(gcopy%m, gcopy%n, alen, gcopy%irn, gcopy%jptr, err)
  __QRM_INFO_CHECK(err, name,'qrm_colamd',9999)

  ! ccolamd return a 0-based permutation (re-argh!)
  cperm = gcopy%jptr(1:graph%n)+1

  call qrm_dealloc(gcopy%jptr)
  
9999 continue
  ! cleanup and return
  call _qrm_spmat_destroy(gcopy, all=.true., info=err2)

  if(present(info)) info = merge(err, err2, err.ne.0)
  return

end subroutine _qrm_do_colamd
