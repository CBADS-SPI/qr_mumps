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
/** @file qrm_mumps.h
 * Header file for the C interface
 *
 * $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
 * $Author: abuttari $
 * $Version: 2.0 $
 * $Revision: 2244 $
 *
 */
/*##############################################################################################*/

#if defined(zprec) || defined(cprec)
#include <complex.h>
#endif

#include <string.h>
#include "qrm_common_c.h"

struct _qrm_spmat_type_c{
  int          *irn, *jcn;
  _qrm_data_c  *val;
  int          m, n, nz;
  int          *cperm_in;
  int          icntl[20];
  double       rcntl[10];
  long int     gstats[10];
  int          h; 
};


int  _qrm_spmat_init_c(struct _qrm_spmat_type_c *qrm_spmat_c);
int  _qrm_spmat_destroy_c(struct _qrm_spmat_type_c *qrm_spmat_c);
void _qrm_readmat_c(char *matfile, struct _qrm_spmat_type_c *qrm_spmat_c);
int  _qrm_analyse_c(struct _qrm_spmat_type_c *qrm_spmat_c, const char transp);
int  _qrm_factorize_c(struct _qrm_spmat_type_c *qrm_spmat_c, const char transp);
int  _qrm_solve_c(struct _qrm_spmat_type_c *qrm_spmat_c, const char transp,
                  _qrm_data_c *b, _qrm_data_c *x, const int nrhs);
int  _qrm_apply_c(struct _qrm_spmat_type_c *qrm_spmat_c, const char transp,
                  _qrm_data_c *b, const int nrhs);
void _qrm_matmul_c(struct _qrm_spmat_type_c *qrm_spmat_c, const char transp,
                   const _qrm_data_c alpha, _qrm_data_c *x, 
                   const _qrm_data_c beta, _qrm_data_c *y, 
                   const int nrhs);
int  _qrm_matnrm_c(struct _qrm_spmat_type_c *qrm_spmat_c, const char ntype, 
                   _qrm_real_c *nrm);
int  _qrm_vecnrm_c(const _qrm_data_c *x, const int n, const int nrhs, 
                   const char ntype, _qrm_real_c *nrm);
int  _qrm_least_squares_c(struct _qrm_spmat_type_c *qrm_spmat_c, _qrm_data_c *b, 
                          _qrm_data_c *x, const int nrhs);
int  _qrm_min_norm_c(struct _qrm_spmat_type_c *qrm_spmat_c, _qrm_data_c *b, 
                          _qrm_data_c *x, const int nrhs);
int  _qrm_residual_norm_c(struct _qrm_spmat_type_c *qrm_spmat_c, _qrm_data_c *b, 
                          _qrm_data_c *x, const int nrhs, _qrm_real_c *nrm);
int  _qrm_residual_orth_c(struct _qrm_spmat_type_c *qrm_spmat_c, _qrm_data_c *r, 
                          const int nrhs, _qrm_real_c *nrm);


int  _qrm_pseti_c(struct _qrm_spmat_type_c *qrm_spmat_c, const char *string, int val);
int  _qrm_pgeti_c(struct _qrm_spmat_type_c *qrm_spmat_c, const char *string, int *val);
int  _qrm_pgetii_c(struct _qrm_spmat_type_c *qrm_spmat_c, const char *string, long long *val);

