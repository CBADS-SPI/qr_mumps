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
!> @file qrm_symbolic.F90
!! This files contains the routine that does the symbolic factorization.
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief This subroutine computes the symbolic QR factorization of a matrix.


!> This routine completely characterizes the structure of fronts in
!! the elimination tree and does a number of other symbolic operations
!! that are essential for the subsequent numerical
!! factorization. Specifically:
!! - it computes the number of rows in each front (note that the
!!   number of columns was already computed by the @link
!!   _qrm_rowcount @endlink routine)
!! - it computes the column indices for each front
!! - it defines a layer in the tree below which, subtrees are treated
!!   sequentially
!! - it estimates the number of flops done in the numerical factorization
!! - it computes the number of nnz in R and H
!! 
!! @param[in] graph This is the adjacency graph of the matrix to be
!!                  factorized in CSC format. On exit, its adata
!!                  member will be modfied. This is the global data
!!                  structure holding all the information computed in
!!                  the analysis phase and needed for the numerical
!!                  factorization. On input adata\%rc, adata\%cperm,
!!                  adata\%parent, adata\%cp_ptr, adata\%nnodes,
!!                  adata\%icperm, adata\%rperm, adata\%child,
!!                  adata\%childptr must be as produced by @link
!!                  qrm_compress_data @endlink. On output, the result
!!                  will be sotred in the following fields:
!!                  - nfrows: the number of rows in each front
!!                  - fcol and fcol_ptr: all the column indices of
!!                    front i will be stored in fcol(fcol_ptr(i):fcol_ptr(i+1)-1)
!!                  - small: if small(i)=1 the it means that front i is the root
!!                    of a subtree that will be treated sequentially during the
!!                    numerical factorization.
!!
subroutine _qrm_symbolic(graph, info)

  use _qrm_spmat_mod
  use qrm_adata_mod
  use qrm_error_mod
  use qrm_sort_mod
  use qrm_common_mod
  use _qrm_utils_mod
  implicit none

  type(_qrm_spmat_type), target   :: graph
  integer, optional               :: info
  
  integer                         :: i, j, f, p, pp, ppp, root, node, roff, ne, np
  integer                         :: first, c, ib, nlz, nth, leaves, totleaves
  integer                         :: m, n, k, cyc, mb, nb, fm, fn, fk, bh
  real(kind(1.d0))                :: rm, rk, rn, smallth
  integer, allocatable            :: col_map(:), mark(:), stair(:), icperm(:)
  integer, pointer                :: porder(:), rc(:), parent(:), fcol(:), fcol_ptr(:)
  type(_qrm_spmat_type)           :: g_csr
  logical                         :: storer, storeh
  type(qrm_adata_type), pointer   :: adata
  integer(kind=8)                 :: hsize, rsize, fhsize, frsize, fflops

  real(kind(1.d0))                :: flop, speedup

  ! error management
  integer                         :: err
  character(len=*), parameter     :: name='qrm_symbolic'

  err = 0

  ! just to simplify
  adata    => graph%adata
  porder   => adata%cperm
  rc       => adata%rc
  parent   => adata%parent
  
  call qrm_alloc(icperm, graph%n, err)
  __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)
  do i=1, graph%n
     icperm(adata%cperm(i)) = i
  end do

  call qrm_alloc(adata%fcol_ptr, adata%nnodes+1, err)
  __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)
  call qrm_alloc(adata%fcol, sum(rc), err)
  __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)

  ! just to simplify
  fcol     => adata%fcol
  fcol_ptr => adata%fcol_ptr

  
  ! first, determine the columns in each front.  

  ! col_map(j)=f means that global column j is a principal variable in
  ! front f
  call qrm_alloc(col_map, graph%n, err)
  __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)
  call qrm_alloc(mark, adata%nnodes, err)
  __QRM_INFO_CHECK(err, name,'qrm_alloc',9999)
  mark = 0

  fcol_ptr(1:2)=1

  ! extract the first iteration in order to initialize fcol_ptr as
  ! well
  f = 1
  do p=adata%cp_ptr(f), adata%cp_ptr(f+1)-1
     j = porder(p)
     col_map(j) = f
  end do

  do f=2, adata%nnodes
     fcol_ptr(f+1) = fcol_ptr(f)+max(rc(f-1),0)
     do p=adata%cp_ptr(f), adata%cp_ptr(f+1)-1
        j = porder(p)
        col_map(j) = f
     end do
  end do
  
#if defined(debug)
  if(p .ne. graph%n+1) then
     __QRM_PRNT_DBG('("Error in symbolic. i .ne. n ",i5,2x,i5)')p, graph%n
  end if
#endif
  
  ! on input the graph is in csc format. we also need a csr
  ! representation in order to determine the row-subtrees of A'*A
  call _qrm_spmat_convert(graph, g_csr, 'csr', .false., info=err)
  __QRM_INFO_CHECK(err, name,'qrm_spmat_convert',9999)

  do i=1,g_csr%nz
     g_csr%jcn(i) = icperm(g_csr%jcn(i))
  end do

  ! sort the coefficients of the matrix along each row. This will
  ! speedup things below
  call _qrm_sort_mat(g_csr, values=.false., info=err)
  __QRM_INFO_CHECK(err, name,'qrm_sort_mat',9999)

  do i=1,g_csr%nz
     g_csr%jcn(i) = adata%cperm(g_csr%jcn(i))
  end do

 
  ! for every front f, for every variable i in f, determine the
  ! coefficients j in (A'*A)(i,1:i) and 
  cyc=0
  do f=1, adata%nnodes
     do p=adata%cp_ptr(f), adata%cp_ptr(f+1)-1
        i = porder(p)

        ! the the coefficient j exists in row i of A'*A iff i and j
        ! are both present in a row of A. So, for each nnz k in column j of A,
        ! every nnz (k,j) in A defines a coefficient (i,j) in A'*A
        fcol(fcol_ptr(f+1)) = i
        fcol_ptr(f+1) = fcol_ptr(f+1)+1

        do pp = graph%jptr(i), graph%jptr(i+1)-1
           k = graph%irn(pp)
           
           ! for every nnz k in column i, go along the corresponding
           ! row and, in the tree, go up until node f (containing i)
           do ppp=g_csr%iptr(k), g_csr%iptr(k+1)-1
              j = g_csr%jcn(ppp)
              if(icperm(j) .ge. icperm(i)) exit

              ! (i,j) is in tril(A'*A). Go up the tree until node f, for
              ! every node met, add column i to the corresponding
              ! front
              node = col_map(j)
              do
                 ! go up the tree
                 if((mark(node) .eq. i) .or. (node .eq. f)) exit

                 fcol(fcol_ptr(node+1)) = i
                 fcol_ptr(node+1) = fcol_ptr(node+1)+1

                 mark(node) = i
                 node = parent(node)
              end do
           end do
        end do
        
     end do
  end do
  if(err.eq.0) call qrm_dealloc(mark,                     err)
  if(err.eq.0) call qrm_dealloc(icperm,                   err)
  if(err.eq.0) call qrm_alloc(adata%nfrows, adata%nnodes, err)
  if(err.eq.0) call qrm_alloc(adata%weight, adata%nnodes, err) 
  if(err.eq.0) call qrm_alloc(adata%asize,  adata%nnodes, err) 
  if(err.eq.0) call qrm_alloc(adata%csize,  adata%nnodes, err) 
  if(err.eq.0) call qrm_alloc(stair,        maxval(rc)+1, err)
  __QRM_INFO_CHECK(err, name,'qrm_(de)alloc',9999)
  adata%weight = 0

  hsize    = 0  
  rsize    = 0  

  storer = graph%icntl(qrm_keeph_).ge.0
  storeh = graph%icntl(qrm_keeph_).ge.1
  
  ! determine structure and weight of all nodes
  do f=1, adata%nnodes

#if defined(debug)
     col_map=0
#endif

     ! under the assumption of postordered nodes, we can simply sweep
     ! the list of fronts

     ! build the col_map for front f. col_map(k)=j means that global
     ! column k is column j inside front f
     do j=1, rc(f)
        k = fcol(fcol_ptr(f)+j-1)
        col_map(k)=j
     end do

     if(f .eq. 1) then
        roff = 1
     else
        roff = adata%stair(f-1)+1
     end if
     
     stair(1:rc(f)) = 0
     ! assemble the rows from the original matrix
     do p=roff, adata%stair(f)
        ! i is a row of the original matrix to be assembled into front f
        i = adata%rperm(p) 
        ! sweep this row and determine its first coefficient
        first = col_map(g_csr%jcn(g_csr%iptr(i)))
        stair(first) = stair(first)+1
     end do
     
     ! assemble the CBs from the children
     do ppp=adata%childptr(f), adata%childptr(f+1)-1
        c = adata%child(ppp)
        
        ! ne is the number of Householder vectors computed on the
        ! child c. np is the number of fully assembled pivots in c
        ne = min(rc(c), adata%nfrows(c))
        np = adata%cp_ptr(c+1)-adata%cp_ptr(c)
        
        ! count in all the rows on the CB of c
        do i=np+1, ne
           j = fcol(fcol_ptr(c)+i-1)
           first = col_map(j)
           stair(first) = stair(first)+1
        end do
     end do

     ! finalize stair
     do i=2, rc(f)
        stair(i) = stair(i)+stair(i-1)
     end do

     if(rc(f) .gt. 0) then
        adata%nfrows(f) = stair(rc(f))
     else
        adata%nfrows(f) = 0
     end if

     ! determine mb, nb, ib and bh
     call qrm_get(graph, 'qrm_mb', mb)
     call qrm_get(graph, 'qrm_nb', nb)
     call qrm_get(graph, 'qrm_ib', ib)
     call qrm_get(graph, 'qrm_bh', bh)
     
     ! At this point it is possible to determine the computational
     ! weight of node f
     ne = min(rc(f),adata%nfrows(f))
     np = adata%cp_ptr(f+1)-adata%cp_ptr(f)

     do i=1, rc(f)
        stair(i) = max(stair(i),min(i,adata%nfrows(f)))
     end do

     call qrm_get_front_mem(adata%nfrows(f), adata%rc(f), mb, nb, ib, &
          & bh, np, adata%asize(f), adata%csize(f), fhsize, frsize, &
          & storer, storeh, stair=stair)

#if defined (flopscount)    
     fflops = qrm_compute_task_flops('geqrt', adata%nfrows(f), adata%rc(f), &
          & -1, -1, stair, 0, 0)
#else
     call qrm_get_front_flops(adata%nfrows(f), adata%rc(f), stair, &
          & mb, nb, ib, fflops)
#endif
     hsize = hsize+fhsize
     rsize = rsize+frsize
     adata%weight(f) = adata%weight(f)+fflops
     p = parent(f)
     if(p.ne.0) adata%weight(p)= adata%weight(p)+ adata%weight(f)
  end do

  graph%gstats(qrm_e_facto_flops_) = adata%weight(adata%nnodes)
  graph%gstats(qrm_e_nnz_r_) = rsize
  graph%gstats(qrm_e_nnz_h_) = hsize
 
#if defined(debug)
  __QRM_PRNT_DBG('("Total estimated number of MFLOPS: ",i10)') adata%weight(adata%nnodes)
#endif


9999 continue
  ! cleanup and return
  call _qrm_spmat_destroy(g_csr, all=.true.)
  call qrm_dealloc(col_map)
  call qrm_dealloc(stair)

  if(err.gt.0) then
     call qrm_dealloc(mark          , err) 
     call qrm_dealloc(icperm        , err) 
     call qrm_dealloc(adata%fcol    , err) 
     call qrm_dealloc(adata%fcol_ptr, err) 
     call qrm_dealloc(adata%nfrows  , err) 
     call qrm_dealloc(adata%asize   , err) 
     call qrm_dealloc(adata%csize   , err) 
  end if

  if(present(info)) info = err
  return
  
end subroutine _qrm_symbolic
