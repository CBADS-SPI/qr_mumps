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


#ifndef __QRM_COMMON_C_H__
#define __QRM_COMMON_C_H__

double qrm_swtime();
int qrm_gseti_c(const char *string, int val);
int qrm_ggeti_c(const char *string, int *val);
int qrm_ggetii_c(const char *string, long long *val);


int  qrm_init_c(int nthreads);
void qrm_finalize_c();

enum icntl{ 
  qrm_ordering_=0,
  qrm_sing_,
  qrm_minamalg_,
  qrm_mb_,
  qrm_nb_,
  qrm_ib_,
  qrm_bh_,
  qrm_keeph_,
  qrm_rhsnb_};

enum rcntl{ 
  qrm_amalgthr_=0,
  qrm_mem_relax_};

enum ords{
  qrm_auto=0,
  qrm_natural_,
  qrm_given_,
  qrm_colamd_,
  qrm_metis_,
  qrm_scotch_};

enum gstats{
  qrm_e_facto_flops_=0,
  qrm_e_nnz_r_,
  qrm_e_nnz_h_,
  qrm_facto_flops_,
  qrm_nnz_r_,
  qrm_nnz_h_,
  qem_e_facto_mempeak_
};

enum yn{
  qrm_no_=0,
  qrm_yes_};

#endif /* __QRM_COMMON_C_H__ */
