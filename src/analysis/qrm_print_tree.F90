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
!> @file qrm_print_tree.F90
!! This file contains various routine for printing trees
!!
!! $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!! $Author: abuttari $
!! $Version: 2.0 $
!! $Revision: 2244 $
!!
!! ##############################################################################################


subroutine qrm_print_nsteps_tree(file, adata, small)

  use qrm_adata_mod
  use qrm_mem_mod
  implicit none

  type(qrm_adata_type) :: adata
  character :: file*(*)
  logical, optional :: small

  integer :: i, p, n, f
  integer, allocatable :: tmp(:)
  logical :: ismall

  if(present(small)) then
     ismall = small
  else
     ismall = .false.
  end if
  
  open(4, file=file, action='write')
  
  write(4,'("graph G {")')
  write(4,'("node [color=black,")')
  write(4,'("fillcolor=white,")')
  write(4,'("shape=circle,")')
  write(4,'("style=filled")')
  write(4,'("];")')

  do i=1, adata%nnodes

     f = adata%torder(i)
     
     if(adata%small(f).lt.0) then
        if(ismall) then
           write(4,'("node",i6.6,"[fillcolor=gray , label="" node:",i6,"\n m:",i6,"\n n:",'//&
                &'i6,"\n np:",i6,"""];")')f,f,adata%nfrows(f),adata%rc(f),&
                &adata%cp_ptr(f+1)-adata%cp_ptr(f)
        end if
     else if(adata%small(f).gt.0) then
        write(4,'("node",i6.6,"[fillcolor=gray , label="" node:",i6,"\n m:",i6,"\n n:",'//&
             &'i6,"\n np:",i6,"""];")')f,f,adata%nfrows(f),adata%rc(f),&
             &adata%cp_ptr(f+1)-adata%cp_ptr(f)
     else if(adata%small(f).eq.0) then
        write(4,'("node",i6.6,"[fillcolor=white, label="" node:",i6,"\n m:",i6,"\n n:",'//&
             &'i6,"\n np:",i6,"""];")')f,f,adata%nfrows(f),adata%rc(f),&
             &adata%cp_ptr(f+1)-adata%cp_ptr(f)
     end if

  end do



  do i=1, adata%nnodes

     f = adata%torder(i)

     p = adata%parent(f)

     if(p.eq.0) cycle
     
     if((adata%small(p).eq.0)) then
        write(4,'("node",i6.6," -- node",i6.6)')p,f
     else
        if(ismall) then
           write(4,'("node",i6.6," -- node",i6.6)')p,f
        end if
     end if

  end do

  write(4,'("}")')

  close(4)


  return

end subroutine qrm_print_nsteps_tree




