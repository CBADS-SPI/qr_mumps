! /* ##############################################################################################
! **
! ** Copyright 2012 CNRS, INPT
! **  
! ** This file is part of qr_mumps.
! **  
! ** qr_mumps is free software: you can redistribute it and/or modify
! ** it under the terms of the GNU Lesser General Public License as 
! ** published by the Free Software Foundation, either version 3 of 
! ** the License, or (at your option) any later version.
! **  
! ** qr_mumps is distributed in the hope that it will be useful,
! ** but WITHOUT ANY WARRANTY; without even the implied warranty of
! ** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! ** GNU Lesser General Public License for more details.
! **  
! ** You can find a copy of the GNU Lesser General Public License
! ** in the qr_mumps/doc directory.
! **
! ** ##############################################################################################*/
!  
!  
! /*##############################################################################################*/
! /** @file qrm_common.h
!  * Common header file
!  *
!  * $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
!  * $Author: abuttari $
!  * $Version: 2.0 $
!  * $Revision: 2244 $
!  *
!  */
! /* ############################################################################################## */

#define __QRM_PRNT_ERR(X)  if(qrm_eunit.gt.0) write(qrm_eunit,X)
#define __QRM_PRNT_MSG(X)  if(qrm_ounit.gt.0) write(qrm_ounit,X)
#define __QRM_PRNT_DBG(X)  if(qrm_dunit.gt.0) write(qrm_dunit,X)




#define __QRM_INFO_CHECK(INFO,NAME,STR,GOTO) if(INFO .ne. 0 ) then; \
 call qrm_error_print(17,NAME,aed=STR,ied=(/INFO/));               \
goto GOTO;\
endif

#define __QRM_INFO_SET(CODE,NAME,GOTO)  info=CODE; \
 call qrm_error_print(CODE,NAME); \
goto GOTO;

#define __QRM_INFO_SET_I(CODE,NAME,IED,GOTO)  info=CODE; \
 call qrm_error_print(CODE,NAME,ied=IED); \
goto GOTO;

#define __QRM_INFO_SET_A(CODE,NAME,AED,GOTO)  info=CODE; \
 call qrm_error_print(CODE,NAME,aed=AED); \
goto GOTO;


#define __QRM_STATUS_SET(CODE,WHERE,GOTO) if(info .ne. 0 ) then;     \
 call qrm_status_set(qrm_dscr,CODE,WHERE);                              \
goto GOTO;\
endif

#define __QRM_SSET_GOTO_I(CODE,WHERE,IED,GOTO) if(info .ne. 0 ) then;  \
 call qrm_status_set(qrm_dscr,CODE,WHERE,ied=IED);                             \
goto GOTO;\
endif
