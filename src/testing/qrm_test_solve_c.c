/* ##############################################################################################
**
** Copyright 2012 CNRS, INPT
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


/**##############################################################################################
* @file qrm_test_solve_c.F90
* This file contains coverage tests for the C interface
*
* $Date: 2016-06-29 14:39:53 +0200 (Wed, 29 Jun 2016) $
* $Author: abuttari $
* $Version: 2.0 $
* $Revision: 2244 $
*
##############################################################################################*/

#include "_qrm_mumps.h"
#include <string.h>
#include <stdio.h>

int _qrm_test_solve_c(struct _qrm_spmat_type_c *qrm_spmat, _qrm_data_c *b, 
                      _qrm_data_c *x, _qrm_data_c *r, _qrm_real_c eps, _qrm_real_c *err){
}



int _qrm_test_methods_c(struct _qrm_spmat_type_c *qrm_spmat, _qrm_data_c *b, 
                        _qrm_data_c *x, _qrm_data_c *r, _qrm_real_c eps, _qrm_real_c *err){

  _qrm_real_c xnrm, bnrm, anrm, rnrm;
  int i, info;

  for(i=0; i<qrm_spmat->m; i++)
    r[i] = b[i];
  info = _qrm_vecnrm(b, qrm_spmat->m, 1, '2', bnrm); 
  if(info > 0) return info;

  if(qrm_spmat->m >= qrm_spmat->n) {
    info  = _qrm_least_squares_c(qrm_spmat, b, x, 1);
  } else {
    info = _qrm_min_norm_c(qrm_spmat, b, x, 1);
  }
  if(info > 0) return info;
  
  info = _qrm_residual_norm_c(qrm_spmat, r, x, 1, &rnrm);
  info = _qrm_residual_orth_c(qrm_spmat, r, 1, &onrm);
  info = _qrm_matnrm(qrm_spmat, 'f', &anrm); 
  info = _qrm_vecnrm(x, qrm_spmat->n, 1, '2', xnrm); 
  if(info > 0) return info;

  if(qrm_spmat->m >= qrm_spmat->n) {
    if(rnrm < eps) {
      *err = rnrm/anrm;
    } else {
      *err = onrm/rnrm;
    }
  } else {
    *err = rnrm/anrm;
  }
  
  if(err < err1)
    err = err1;

  return info;

}
