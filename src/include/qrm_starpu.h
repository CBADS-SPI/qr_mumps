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


#ifndef __QRM_STARPU_H__
#define __QRM_STARPU_H__

enum qrm_starpu_op {
  activate,
  panel,
  update,
  assemble,
  clean,
  callback,
  idle,
  update_gpu
};

/** 2D kernels **/

struct pnl_params {
  unsigned nb;
  unsigned pnl;
  unsigned parts;
  /* inner blocking for tiled factorization*/
  unsigned ib;
};

struct upd_params {
  unsigned nb;
  unsigned pnl;
  unsigned pnlUp;
  /* inner blocking for tiled factorization*/
  unsigned ib;
};

#define MIN(a,b)	((a)<(b)?(a):(b))
#define MAX(a,b)	((a)<(b)?(b):(a))

#if defined(usegpu)
#define MAX_STREAM 1
#endif

/* performance model */

double expected_time_nl_regression_outer_upd_cpu(double size);
double expected_time_nl_regression_outer_upd_cuda(double size);

#endif /* __QRM_STARPU_H__ */
