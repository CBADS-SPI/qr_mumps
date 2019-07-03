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


#include <starpu.h>
#include <starpu_cuda.h>
#include <starpu_profiling.h>
#include <limits.h>

int starpu_f_init_c(int ncpus) {

  int info;
  struct starpu_conf *conf = malloc(sizeof(struct starpu_conf));
  
  starpu_conf_init(conf);

  if(ncpus > 0)
    conf->ncpus = ncpus;  

  /* conf->sched_policy_name = "eager"; */

  info = starpu_init(conf);

  free(conf);
  return info;
}

void starpu_f_get_buffer(void *buffers[], int num, void **A, int *m, int *n, int *lda) {

  *A   = (void *)STARPU_MATRIX_GET_PTR(buffers[num]);
  *m   = (int)STARPU_MATRIX_GET_NX(buffers[num]);
  *n   = (int)STARPU_MATRIX_GET_NY(buffers[num]);
  *lda = (int)STARPU_MATRIX_GET_LD(buffers[num]);

  return;

}


void starpu_f_data_acquire_read(starpu_data_handle_t data) {

  starpu_data_acquire(data, STARPU_R);
  
  return;

}


int starpu_f_sched_ctx_create_c(int *workers, int nworkers, const char *name){
  unsigned ctx;

  int i;
  
  ctx = starpu_sched_ctx_create(workers, nworkers, name, STARPU_SCHED_CTX_POLICY_NAME, NULL, 0);
  
  return (int)ctx;
}


void starpu_f_sched_ctx_set_context_c(int *ctx){

  starpu_sched_ctx_set_context((unsigned*)ctx);

}


void starpu_f_sched_ctx_display_workers_c(int ctx){

  starpu_sched_ctx_display_workers((unsigned)ctx, stderr);

}

void starpu_f_sched_ctx_delete_c(int ctx){

  starpu_sched_ctx_delete((unsigned)ctx);

}

void starpu_f_task_wait_for_all_in_ctx_c(int ctx){

  starpu_task_wait_for_all_in_ctx((unsigned)ctx);

}
