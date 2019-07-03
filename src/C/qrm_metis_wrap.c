/* ##############################################################################################
**
** Copyright 2012-2016 CNRS, INPT
** Copyright 2013-2015 UPS
**  
** This file is part of qr_mumps.
**  
** qr_mumps is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License as 
** published by the Free Software Foundation, either version 3 of 
** the License, or (at your option) any later version.
**  
** qr_mumps is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**  
** You can find a copy of the GNU Lesser General Public License
** in the qr_mumps/doc directory.
**
** ##############################################################################################*/


/*##############################################################################################*/
/** @file qrm_metis_wrap.c
 * FIXME: add comments
 *
 * $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
 * $Author: abuttari $
 * $Version 0.0.1$
 * $Revision: 2244 $
 *
 **/
/*##############################################################################################*/


#if defined(have_metis)
#include "metis.h"
#endif
#include <stdio.h>

void qrm_metis(int *n, int *iptr, int *jcn, int *cperm, int *iperm, int *info){
#if defined(have_metis)


#if defined(METIS_VER_MAJOR)
  int options [METIS_NOPTIONS] ;	
  int err;
  
  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_NUMBERING]=1;

  err = METIS_NodeND(n, iptr, jcn, NULL, options, cperm, iperm);

  if(err != METIS_OK) *info=20;
  
#else
  int options [8] ;	
  int numflag;

  options[0]=0;

  numflag=1;
  METIS_NodeND(n, iptr, jcn, &numflag, &options[0], cperm, iperm);

  *info = 0;
#endif

  return;
#endif
}







