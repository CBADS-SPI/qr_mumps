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


#include "_qrm_c.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

void readmat(char *fname, struct _qrm_spmat_type_c *qrm_mat);

int main(){
  struct _qrm_spmat_type_c qrm_mat, rmat;
  int i, kh;
  char transp;
  _qrm_data_c *b, *x, *r;
  double t1, ta, tf, ts;
  _qrm_real_c rnrm, onrm, anrm, bnrm, xnrm;

  qrm_init_c(0);
  
  printf("==================================\n");
  _qrm_spmat_init_c(&qrm_mat);

  printf("==================================\n");

  readmat("image_interp.mtx", &qrm_mat);

  qrm_gseti_c("qrm_dunit", -1);


  /* qrm_mat.icntl[QRM_ORDERING_] = QRM_SCOTCH_; */
  _qrm_pseti_c(&qrm_mat, "qrm_ordering",qrm_metis_);

  /* qrm_mat.icntl[QRM_KEEPH_]    = QRM_YES_; */
  _qrm_pseti_c(&qrm_mat, "qrm_keeph",qrm_yes_);


  if(qrm_mat.m < qrm_mat.n){
    transp = _qrm_transp;
  } else {
    transp = 'n';
  }

  printf("\nStarting the analysis\n");
  t1 = qrm_swtime();
  _qrm_analyse_c(&qrm_mat, transp);
  ta = qrm_swtime()-t1;

  printf("Starting the factorization\n");
  t1 = qrm_swtime();
  _qrm_factorize_c(&qrm_mat, transp);
  tf = qrm_swtime()-t1;

  _qrm_pgeti_c(&qrm_mat, "qrm_keeph", &kh);

  if(kh == qrm_yes_ ){
    x = (_qrm_data_c *)malloc(qrm_mat.n*sizeof(_qrm_data_c));
    b = (_qrm_data_c *)malloc(qrm_mat.m*sizeof(_qrm_data_c));
    r = (_qrm_data_c *)malloc(qrm_mat.m*sizeof(_qrm_data_c));
    
    for(i=0; i<qrm_mat.m; i++){
      b[i] = _qrm_one_c;
      r[i] = _qrm_one_c;
    }

    for(i=0; i<qrm_mat.n; i++){
      x[i] = _qrm_zero_c;
    }

    printf("Starting the solve\n");
    t1 = qrm_swtime();
    if(transp == 'n'){
      _qrm_apply_c(&qrm_mat, _qrm_transp, b, 1);
      _qrm_solve_c(&qrm_mat, 'n', b, x, 1);
    } else if (transp == _qrm_transp) {
      _qrm_solve_c(&qrm_mat, _qrm_transp, b, x, 1);
      _qrm_apply_c(&qrm_mat, 'n', x, 1);
    }
    ts = qrm_swtime()-t1;

    printf("Error checking\n");

    _qrm_residual_norm_c(&qrm_mat, r, x, 1, &rnrm);
    _qrm_residual_orth_c(&qrm_mat, r, 1, &onrm);
    _qrm_vecnrm_c(x, qrm_mat.n, 1, '2', &xnrm);
    _qrm_vecnrm_c(b, qrm_mat.m, 1, '2', &bnrm);
    _qrm_matnrm_c(&qrm_mat, 'f', &anrm);

    printf("||A||          =  %10.5e\n",anrm);
    printf("||b||          =  %10.5e\n",bnrm);
    printf("||x||          =  %10.5e\n",xnrm);
    printf("||r||/||A||    =  %10.5e\n",rnrm);
    printf("||A^tr||/||r|| =  %10.5e\n",onrm);

    
    free(x); free(b); free(r);
  }

  printf("\n\nTime to do the analysis      : %e\n",ta);
  printf("Time to do the factorization : %e\n",tf);
  printf("Time to do the solve         : %e\n",ts);
  printf("Nonzeros in R                : %ld\n",qrm_mat.gstats[qrm_nnz_r_]);
  printf("Nonzeros in H                : %ld\n",qrm_mat.gstats[qrm_nnz_h_]);
  printf("Total flops at facto         : %ld\n",qrm_mat.gstats[qrm_e_facto_flops_]);
  
  _qrm_spmat_destroy_c(&qrm_mat);
  free(qrm_mat.irn); free(qrm_mat.jcn); free(qrm_mat.val);

  qrm_finalize_c();
  return 1;

}


#define MM_MAX_LINE_LENGTH 1025
#define MatrixMarketBanner "%%MatrixMarket"
#define MM_MAX_TOKEN_LENGTH 64

void readmat(char *fname, struct _qrm_spmat_type_c *qrm_mat){
  /* Largely copied from mmio.c from Matrix Market */

  FILE * file;
  file = fopen( fname, "r");
  char line[MM_MAX_LINE_LENGTH];
  char banner[MM_MAX_TOKEN_LENGTH];
  char mtx[MM_MAX_TOKEN_LENGTH]; 
  char crd[MM_MAX_TOKEN_LENGTH];
  char data_type[MM_MAX_TOKEN_LENGTH];
  char storage_scheme[MM_MAX_TOKEN_LENGTH];
  char *p;
  int i;
  _qrm_real_c re, im;
  
  if (fgets(line, MM_MAX_LINE_LENGTH, file) == NULL) {
    printf("Error!\n");
    return;
  }

  if (sscanf(line, "%s %s %s %s %s", banner, mtx, crd, data_type, 
	     storage_scheme) != 5) {
    printf("Error!\n");
    return;
  }

  for (p=mtx; *p!='\0'; *p=tolower(*p),p++);  /* convert to lower case */
  for (p=crd; *p!='\0'; *p=tolower(*p),p++);  
  for (p=data_type; *p!='\0'; *p=tolower(*p),p++);
  for (p=storage_scheme; *p!='\0'; *p=tolower(*p),p++);
  
  /* check for banner */
  if (strncmp(banner, MatrixMarketBanner, strlen(MatrixMarketBanner)) != 0){
    printf("Error!\n");
    return;
  }

  
  do {
    if (fgets(line,MM_MAX_LINE_LENGTH,file) == NULL) {
      printf("Error!\n");
      return;
    }
  } while (line[0] == '%');
  
  sscanf(line, "%d %d %d", &qrm_mat->m, &qrm_mat->n, &qrm_mat->nz);

  printf("Reading Matrix %s -- m=%d  n=%d  nnz=%d\n", fname, qrm_mat->m, qrm_mat->n, qrm_mat->nz);
  
  qrm_mat->irn = (int *)malloc(qrm_mat->nz*sizeof(int));
  qrm_mat->jcn = (int *)malloc(qrm_mat->nz*sizeof(int));
  qrm_mat->val = (_qrm_data_c *)malloc(qrm_mat->nz*sizeof(_qrm_data_c));

  
  if (strncmp(data_type, "pattern", 7) == 0){
    for (i=0; i<qrm_mat->nz; i++) {
      fscanf(file, "%d %d\n", &qrm_mat->irn[i], &qrm_mat->jcn[i]);
      qrm_mat->val[i] = _qrm_one_c;
    }
  } else {
    for (i=0; i<qrm_mat->nz; i++) {
#if defined(dprec)
      fscanf(file, "%d %d %lg\n", &qrm_mat->irn[i], &qrm_mat->jcn[i], &qrm_mat->val[i]);
#elif defined(sprec)
      fscanf(file, "%d %d %g\n", &qrm_mat->irn[i], &qrm_mat->jcn[i], &qrm_mat->val[i]);
#elif defined(zprec)
      fscanf(file, "%d %d %lg %lg\n", &qrm_mat->irn[i], &qrm_mat->jcn[i], &re, &im);
      qrm_mat->val[i] = re+_Complex_I*im;
#elif defined(cprec)
      fscanf(file, "%d %d %g %g\n", &qrm_mat->irn[i], &qrm_mat->jcn[i], &re, &im);
      qrm_mat->val[i] = re+_Complex_I*im;
#endif
    }
  }


}
