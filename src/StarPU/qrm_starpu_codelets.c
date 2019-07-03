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

#include "qrm_starpu.h"
/* #include "qrm_starpu_prof.h" */

#include <starpu.h>
/* #include <starpu_profiling.h> */
#include <stdio.h>
#include <limits.h>


 /* These two types have to be interoperable. Always make sure they
    are aligned with the F90 counterparts in qrm_pthreads_mod.F90 and
    qrm_dscr_mod.F90 */

typedef struct {
  void *m;
} qrm_pthread_mutex;

typedef struct {
  int err_status;
  qrm_pthread_mutex mutex;
  int ctx;
} qrm_dscr_type;


/*******************************************************************************/
/* do_subtree task */
/*******************************************************************************/
void _qrm_starpu_do_subtree_cpu_func(void *buffers[], void *cl_arg);

struct starpu_codelet cl_do_subtree = {
  .where = STARPU_CPU,
  .cpu_funcs = {_qrm_starpu_do_subtree_cpu_func, NULL},
  .nbuffers = 1,
  .modes = {STARPU_RW},
  .name = "Do_subtree"
};

/* insert task to process subtree */
void _qrm_insert_do_subtree_task_c(qrm_dscr_type *qrm_dscr, void *qrm_mat,
                                   int rootnum,
                                   starpu_data_handle_t sym_handle,
                                   double flops, int prio) {

  int ret;

  /* printf("insert_do_subtree_task_c, rootnum: %d %d\n", rootnum, (int)sym_handle); */

  ret = starpu_task_insert(&cl_do_subtree,
                           STARPU_VALUE,     &qrm_dscr, sizeof(void *),
                           STARPU_VALUE,     &qrm_mat, sizeof(void *),
                           STARPU_VALUE,     &rootnum, sizeof(int),
                           STARPU_RW,        sym_handle,
                           STARPU_PRIORITY,  prio,
                           STARPU_SCHED_CTX, (unsigned)qrm_dscr->ctx,
                           0);
  STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_insert");

  return;
}





/*******************************************************************************/
/* init task */
/*******************************************************************************/
void _qrm_starpu_init_cpu_func(void *buffers[], void *cl_arg);


struct starpu_codelet cl_init = {
  .where     = STARPU_CPU,
  .cpu_funcs = {_qrm_starpu_init_cpu_func, NULL},
  .nbuffers  = STARPU_VARIABLE_NBUFFERS,
  .name      = "INIT"
};


void _qrm_insert_init_task_c(qrm_dscr_type *qrm_dscr, void *qrm_mat,
                             int fnum,
                             starpu_data_handle_t front_handle,
                             starpu_data_handle_t *kids_handles,
                             int nc, int prio) {

  int ret;
  int i;
  struct starpu_data_descr *descrs;
  descrs = malloc((nc+1) * sizeof(struct starpu_data_descr));

  descrs[0].handle = front_handle; descrs[0].mode = STARPU_RW;
  for(i=0; i<nc; i++){
    descrs[i+1].handle = kids_handles[i];  descrs[i+1].mode = STARPU_R;
  }

  ret = starpu_task_insert(&cl_init,
                           STARPU_VALUE,           &qrm_dscr, sizeof(void *),
                           STARPU_VALUE,           &qrm_mat, sizeof(void *),
                           STARPU_VALUE,           &fnum, sizeof(int),
                           STARPU_DATA_MODE_ARRAY, descrs,   nc+1,
                           STARPU_PRIORITY,        prio,
                           STARPU_SCHED_CTX,       (unsigned)qrm_dscr->ctx,
                           0);

  STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_insert");

  free(descrs);

  return;
}



/*******************************************************************************/
/* _geqrt_sym task */
/*******************************************************************************/

void _qrm_starpu_geqrt_cpu_func(void *buffers[], void *cl_arg);

struct starpu_codelet cl_geqrt = {
  .where     = STARPU_CPU,
  .cpu_funcs = {_qrm_starpu_geqrt_cpu_func, NULL},
  .nbuffers  = STARPU_VARIABLE_NBUFFERS,
  .name      = "GEQRT"
};


void _qrm_insert_geqrt_task_c(qrm_dscr_type *qrm_dscr,
                              starpu_data_handle_t a_handle,
                              starpu_data_handle_t t_handle,
                              starpu_data_handle_t work_handle,
                              starpu_data_handle_t sym_handle,
                              int bk, int bi,
                              int mb, int nb, int ib,
                              int* stair, double flops, int prio) {

  int ret;

  /* printf("> insert geqrt  %d %d\n", bk, bi); */
  if(sym_handle){
    ret = starpu_task_insert(&cl_geqrt,
                             STARPU_VALUE,     &bk,         sizeof(int),
                             STARPU_VALUE,     &bi,         sizeof(int),
                             STARPU_VALUE,     &mb,         sizeof(int),
                             STARPU_VALUE,     &nb,         sizeof(int),
                             STARPU_VALUE,     &ib,         sizeof(int),
                             STARPU_VALUE,     &stair,      sizeof(int *),
                             STARPU_RW,        a_handle,
                             STARPU_RW,        t_handle,
                             STARPU_SCRATCH,   work_handle,
                             STARPU_R,         sym_handle,
                             STARPU_PRIORITY,  prio,
                             STARPU_SCHED_CTX, (unsigned)qrm_dscr->ctx,
                             0);
  } else {
    ret = starpu_task_insert(&cl_geqrt,
                             STARPU_VALUE,     &bk,         sizeof(int),
                             STARPU_VALUE,     &bi,         sizeof(int),
                             STARPU_VALUE,     &mb,         sizeof(int),
                             STARPU_VALUE,     &nb,         sizeof(int),
                             STARPU_VALUE,     &ib,         sizeof(int),
                             STARPU_VALUE,     &stair,      sizeof(int *),
                             STARPU_RW,        a_handle,
                             STARPU_RW,        t_handle,
                             STARPU_SCRATCH,   work_handle,
                             STARPU_PRIORITY,  prio,
                             STARPU_SCHED_CTX, (unsigned)qrm_dscr->ctx,
                             0);
  }
  STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_insert");

  return;
}





/*******************************************************************************/
/* gemqrt task */
/*******************************************************************************/

void _qrm_starpu_gemqrt_cpu_func(void *buffers[], void *cl_arg);

struct starpu_codelet cl_gemqrt = {
  .where     = STARPU_CPU,
  .cpu_funcs = {_qrm_starpu_gemqrt_cpu_func, NULL},
  .nbuffers  = STARPU_VARIABLE_NBUFFERS,
  .name      = "GEMQRT",
};


void _qrm_insert_gemqrt_task_c(qrm_dscr_type *qrm_dscr,
                               starpu_data_handle_t t_handle,
                               starpu_data_handle_t c_handle,
                               starpu_data_handle_t work_handle,
                               starpu_data_handle_t sym_handle,
                               int bk, int bi, int bj,
                               int mb, int nb, int ib,
                               int* stair, double flops, int prio) {

  int prio_level, ret;

  /* printf("> insert update, fnum : %d, %d, %d, %d\n", bk, bi, bj, (int)sym_handle); */

  if(sym_handle){
    ret = starpu_task_insert(&cl_gemqrt,
                             STARPU_VALUE,     &bk, sizeof(int),
                             STARPU_VALUE,     &bi, sizeof(int),
                             STARPU_VALUE,     &bj, sizeof(int),
                             STARPU_VALUE,     &mb, sizeof(int),
                             STARPU_VALUE,     &nb, sizeof(int),
                             STARPU_VALUE,     &ib, sizeof(int),
                             STARPU_VALUE,     &stair, sizeof(int *),
                             STARPU_R, 	       t_handle,
                             STARPU_RW,	       c_handle,
                             STARPU_SCRATCH,   work_handle,
                             STARPU_R,         sym_handle,
                             STARPU_PRIORITY,  prio,
                             STARPU_SCHED_CTX, (unsigned)qrm_dscr->ctx,
                             0);
  } else {
    ret = starpu_task_insert(&cl_gemqrt,
                             STARPU_VALUE,     &bk, sizeof(int),
                             STARPU_VALUE,     &bi, sizeof(int),
                             STARPU_VALUE,     &bj, sizeof(int),
                             STARPU_VALUE,     &mb, sizeof(int),
                             STARPU_VALUE,     &nb, sizeof(int),
                             STARPU_VALUE,     &ib, sizeof(int),
                             STARPU_VALUE,     &stair, sizeof(int *),
                             STARPU_R, 	       t_handle,
                             STARPU_RW,	       c_handle,
                             STARPU_SCRATCH,   work_handle,
                             STARPU_PRIORITY,  prio,
                             STARPU_SCHED_CTX, (unsigned)qrm_dscr->ctx,
                             0);
  }

  STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_insert");

  return;
}





/*******************************************************************************/
/* tpqrt task */
/*******************************************************************************/

void _qrm_starpu_tpqrt_cpu_func(void *buffers[], void *cl_arg);

struct starpu_codelet cl_tpqrt = {
  .where     = STARPU_CPU,
  .cpu_funcs = {_qrm_starpu_tpqrt_cpu_func, NULL},
  .nbuffers  = STARPU_VARIABLE_NBUFFERS,
  .name      = "TPQRT",
};

void _qrm_insert_tpqrt_task_c(qrm_dscr_type *qrm_dscr,
                              starpu_data_handle_t a_handle,
			      starpu_data_handle_t b_handle,
			      starpu_data_handle_t t_handle,
                              starpu_data_handle_t work_handle,
                              starpu_data_handle_t sym_handle,
                              int bk, int bi, int bl,
                              int mb, int nb, int ib,
                              char ts, int* stair,
                              double flops, int prio) {

  int ret;

  /* printf("> insert update, fnum : %d, %d->%d\n", fnum, pnl, upd); */

  if(sym_handle){
    ret = starpu_task_insert(&cl_tpqrt,
                             STARPU_VALUE,     &bk,       sizeof(int),
                             STARPU_VALUE,     &bi,       sizeof(int),
                             STARPU_VALUE,     &bl,       sizeof(int),
                             STARPU_VALUE,     &mb,       sizeof(int),
                             STARPU_VALUE,     &nb,       sizeof(int),
                             STARPU_VALUE,     &ib,       sizeof(int),
                             STARPU_VALUE,     &ts,       sizeof(char),
                             STARPU_VALUE,     &stair,    sizeof(int *),
                             STARPU_RW,        a_handle,
                             STARPU_RW,	       b_handle,
                             STARPU_RW,	       t_handle,
                             STARPU_SCRATCH,   work_handle,
                             STARPU_R,         sym_handle,
                             STARPU_PRIORITY,  prio,
                             STARPU_SCHED_CTX, (unsigned)qrm_dscr->ctx,
                             0);
  } else {
    ret = starpu_task_insert(&cl_tpqrt,
                             STARPU_VALUE,     &bk,       sizeof(int),
                             STARPU_VALUE,     &bi,       sizeof(int),
                             STARPU_VALUE,     &bl,       sizeof(int),
                             STARPU_VALUE,     &mb,       sizeof(int),
                             STARPU_VALUE,     &nb,       sizeof(int),
                             STARPU_VALUE,     &ib,       sizeof(int),
                             STARPU_VALUE,     &ts,       sizeof(char),
                             STARPU_VALUE,     &stair,    sizeof(int *),
                             STARPU_RW,        a_handle,
                             STARPU_RW,	       b_handle,
                             STARPU_RW,	       t_handle,
                             STARPU_SCRATCH,   work_handle,
                             STARPU_PRIORITY,  prio,
                             STARPU_SCHED_CTX, (unsigned)qrm_dscr->ctx,
                             0);
  }
  STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_insert");

  return;
}





/*******************************************************************************/
/* tpmqrt task */
/*******************************************************************************/

void _qrm_starpu_tpmqrt_cpu_func(void *buffers[], void *cl_arg);

struct starpu_codelet cl_tpmqrt = {
  .where     = STARPU_CPU,
  .cpu_funcs = {_qrm_starpu_tpmqrt_cpu_func, NULL},
  .nbuffers  = STARPU_VARIABLE_NBUFFERS,
  .name      = "TPMQRT",
};


void _qrm_insert_tpmqrt_task_c(qrm_dscr_type *qrm_dscr,
                               starpu_data_handle_t v_handle,
			       starpu_data_handle_t t_handle,
			       starpu_data_handle_t a_handle,
			       starpu_data_handle_t b_handle,
			       starpu_data_handle_t work_handle,
                               starpu_data_handle_t sym_handle,
                               int bk, int bi, int bl, int bj,
                               int mb, int nb, int ib,
                               char ts, int* stair,
                               double flops, int prio) {

  int ret;

  /* printf("> insert update, fnum : %d, %d->%d\n", fnum, pnl, upd); */

  if(sym_handle){
    ret = starpu_task_insert(&cl_tpmqrt,
                             STARPU_VALUE,     &bk,       sizeof(int),
                             STARPU_VALUE,     &bi,       sizeof(int),
                             STARPU_VALUE,     &bl,       sizeof(int),
                             STARPU_VALUE,     &bj,       sizeof(int),
                             STARPU_VALUE,     &mb,       sizeof(int),
                             STARPU_VALUE,     &nb,       sizeof(int),
                             STARPU_VALUE,     &ib,       sizeof(int),
                             STARPU_VALUE,     &ts,       sizeof(char),
                             STARPU_VALUE,     &stair,    sizeof(int *),
                             STARPU_R,         v_handle,
                             STARPU_R,	       t_handle,
                             STARPU_RW,	       a_handle,
                             STARPU_RW,	       b_handle,
                             STARPU_SCRATCH,   work_handle,
                             STARPU_R,         sym_handle,
                             STARPU_PRIORITY,  prio,
                             STARPU_SCHED_CTX, (unsigned)qrm_dscr->ctx,
                             0);
  } else {
    ret = starpu_task_insert(&cl_tpmqrt,
                             STARPU_VALUE,     &bk,       sizeof(int),
                             STARPU_VALUE,     &bi,       sizeof(int),
                             STARPU_VALUE,     &bl,       sizeof(int),
                             STARPU_VALUE,     &bj,       sizeof(int),
                             STARPU_VALUE,     &mb,       sizeof(int),
                             STARPU_VALUE,     &nb,       sizeof(int),
                             STARPU_VALUE,     &ib,       sizeof(int),
                             STARPU_VALUE,     &ts,       sizeof(char),
                             STARPU_VALUE,     &stair,    sizeof(int *),
                             STARPU_R,         v_handle,
                             STARPU_R,	       t_handle,
                             STARPU_RW,	       a_handle,
                             STARPU_RW,	       b_handle,
                             STARPU_SCRATCH,   work_handle,
                             STARPU_PRIORITY,  prio,
                             STARPU_SCHED_CTX, (unsigned)qrm_dscr->ctx,
                             0);
  }



  STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_insert");

  return;
}





/*******************************************************************************/
/* assembly task */
/*******************************************************************************/

void _qrm_starpu_assemble_cpu_func(void *buffers[], void *cl_arg);


struct starpu_codelet cl_asm = {
  .where     = STARPU_CPU,
  .cpu_funcs = {_qrm_starpu_assemble_cpu_func, NULL},
  .nbuffers  = STARPU_VARIABLE_NBUFFERS,
  .name      = "ASM"
};

/* submit an assembly task
 **/
void _qrm_insert_asm_blk_task_c(qrm_dscr_type *qrm_dscr,
                                void *qrm_mat, int fnum,
                                starpu_data_handle_t block_handle,
                                starpu_data_handle_t sym_handle,
                                starpu_data_handle_t *father_handles,
                                starpu_data_handle_t father_sym_handle,
                                int nh, int br, int bc, int prio) {

  int ret;
  int i;
  struct starpu_data_descr *descrs;

  /* printf("> insert assemble block, fnum: %d, col: %d, fcol: %d\n", fnum, col, fcol); */
  /* printf("> insert assemble block, fathernum: %d\n", fathernum); */

  descrs = malloc((nh+3) * sizeof(struct starpu_data_descr));

  descrs[0].handle = sym_handle;	descrs[0].mode = STARPU_R;
  descrs[1].handle = father_sym_handle; descrs[1].mode = STARPU_R;
  descrs[2].handle = block_handle;	descrs[2].mode = STARPU_R;

  for(i=0; i<nh; i++){
    /* printf("insrt asm %p \n",father_handles[i]); */
    descrs[i+3].handle = father_handles[i];  descrs[i+3].mode = STARPU_RW | STARPU_COMMUTE;
  }

  ret = starpu_task_insert(&cl_asm,
                           STARPU_VALUE,           &qrm_mat, sizeof(void *),
                           STARPU_VALUE,           &fnum,    sizeof(int),
                           STARPU_VALUE,           &br,      sizeof(int),
                           STARPU_VALUE,           &bc,      sizeof(int),
                           STARPU_DATA_MODE_ARRAY, descrs,   nh+3,
                           STARPU_PRIORITY,        prio,
                           STARPU_SCHED_CTX,       (unsigned)qrm_dscr->ctx,
                           0);
  STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_insert");

  /* unsure whether the following free is legal or not */
  free(descrs);

  return;

}


/*******************************************************************************/
/* clean task */
/*******************************************************************************/

void _qrm_starpu_clean_cpu_func(void *buffers[], void *cl_arg);

struct starpu_codelet cl_clean = {
  .where     = STARPU_CPU,
  .cpu_funcs = {_qrm_starpu_clean_cpu_func, NULL},
  .nbuffers  = 1,
  .modes     = {STARPU_RW},
  .name      = "CLEAN",
};


/* submit cleaning task
**/
void _qrm_insert_clean_task_c(qrm_dscr_type *qrm_dscr,
                              void *qrm_mat,
                              int fnum,
                              starpu_data_handle_t sym_handle,
                              int prio){

  int ret;
  /* printf("clean tsk %d\n",(int)sym_handle); */
  ret = starpu_task_insert(&cl_clean,
                           STARPU_VALUE,     &qrm_dscr, sizeof(void *),
                           STARPU_VALUE,     &qrm_mat, sizeof(void *),
                           STARPU_VALUE,     &fnum, sizeof(int),
                           STARPU_RW,        sym_handle,
                           STARPU_PRIORITY,  prio,
                           STARPU_SCHED_CTX, (unsigned)qrm_dscr->ctx,
                           0);
  STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_insert");

  return;

}

/*******************************************************************************/
/* apply task */
/*******************************************************************************/

void _qrm_starpu_apply_cpu_func(void *buffers[], void *cl_arg);


struct starpu_codelet cl_apply = {
  .where     = STARPU_CPU,
  .cpu_funcs = {_qrm_starpu_apply_cpu_func, NULL},
  .nbuffers  = STARPU_VARIABLE_NBUFFERS,
  .name      = "Apply"
};

/* submit an apply task
**/
void _qrm_insert_apply_task_c(qrm_dscr_type *qrm_dscr, void *qrm_mat,
                              void *b,
                              starpu_data_handle_t front_handle,
                              starpu_data_handle_t node_handle,
                              starpu_data_handle_t *kids_handles,
                              starpu_data_handle_t work_handle,
                              char transp, int inode, int nc,
                              int prio) {

  int ret;
  int i;
  struct starpu_data_descr *descrs;

  /* printf("> insert apply, inode:%d  nc:%d  transp:%c\n", inode, nc, transp); */

  descrs = malloc((nc+3) * sizeof(struct starpu_data_descr));

  descrs[0].handle = work_handle ; descrs[0].mode = STARPU_SCRATCH;
  descrs[1].handle = front_handle; descrs[1].mode = STARPU_R;
  if(transp==_qrm_transp){
    descrs[2].handle = node_handle; descrs[2].mode = STARPU_RW;
    for(i=0; i<nc; i++){
      descrs[i+3].handle = kids_handles[i];  descrs[i+3].mode = STARPU_R;
    }
  } else {
    descrs[2].handle = node_handle; descrs[2].mode = STARPU_RW;
    for(i=0; i<nc; i++){
      descrs[i+3].handle = kids_handles[i];  descrs[i+3].mode = STARPU_RW;
    }
  }

  ret = starpu_task_insert(&cl_apply,
                           STARPU_VALUE,           &qrm_dscr, sizeof(void *),
                           STARPU_VALUE,           &qrm_mat,  sizeof(void *),
                           STARPU_VALUE,           &b,        sizeof(void *),
                           STARPU_VALUE,           &inode,    sizeof(int),
                           STARPU_VALUE,           &transp,   sizeof(char),
                           STARPU_DATA_MODE_ARRAY, descrs,    nc+3,
                           STARPU_PRIORITY,        prio,
                           STARPU_SCHED_CTX,       (unsigned)qrm_dscr->ctx,
                           0);
  STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_insert");

  /* unsure whether the following free is legal or not */
  free(descrs);

  return;

}

/*******************************************************************************/
/* solve task */
/*******************************************************************************/

void _qrm_starpu_solve_cpu_func(void *buffers[], void *cl_arg);


struct starpu_codelet cl_solve = {
  .where     = STARPU_CPU,
  .cpu_funcs = {_qrm_starpu_solve_cpu_func, NULL},
  .nbuffers  = STARPU_VARIABLE_NBUFFERS,
  .name      = "Solve"
};

/* submit an solve task
**/
void _qrm_insert_solve_task_c(qrm_dscr_type *qrm_dscr, void *qrm_mat,
                              void *b, void *x,
                              starpu_data_handle_t front_handle,
                              starpu_data_handle_t x_node_handle,
                              starpu_data_handle_t b_node_handle,
                              starpu_data_handle_t *kids_handles,
                              char transp, int inode, int nc,
                              int prio) {

  int ret;
  int i;
  struct starpu_data_descr *descrs;

  /* printf("> insert solve, inode:%d  nc:%d  transp:%c\n", inode, nc, transp); */

  descrs = malloc((2*nc+3) * sizeof(struct starpu_data_descr));

  descrs[0].handle = front_handle; descrs[0].mode = STARPU_R;
  descrs[1].handle = x_node_handle; descrs[1].mode = STARPU_RW;
  descrs[2].handle = b_node_handle; descrs[2].mode = STARPU_R;
  if(transp==_qrm_transp){
    for(i=0; i<nc; i++){
      descrs[i+3].handle = kids_handles[i];  descrs[i+3].mode = STARPU_R;
    }
  } else {
    for(i=0; i<nc; i++){
      descrs[i+3].handle = kids_handles[i];  descrs[i+3].mode = STARPU_RW;
    }
  }
  for(i=0; i<nc; i++){
    descrs[nc+i+3].handle = kids_handles[nc+i];  descrs[nc+i+3].mode = STARPU_R;
  }

  ret = starpu_task_insert(&cl_solve,
                           STARPU_VALUE,           &qrm_dscr, sizeof(void *),
                           STARPU_VALUE,           &qrm_mat,  sizeof(void *),
                           STARPU_VALUE,           &b,        sizeof(void *),
                           STARPU_VALUE,           &x,        sizeof(void *),
                           STARPU_VALUE,           &inode,    sizeof(int),
                           STARPU_VALUE,           &transp,   sizeof(char),
                           STARPU_DATA_MODE_ARRAY, descrs,    2*nc+3,
                           STARPU_PRIORITY,        prio,
                           STARPU_SCHED_CTX,       (unsigned)qrm_dscr->ctx,
                           0);
  STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_insert");

  /* unsure whether the following free is legal or not */
  free(descrs);

  return;

}



/*******************************************************************************/
/* solve task */
/*******************************************************************************/

void _qrm_starpu_analyse_cpu_func(void *buffers[], void *cl_arg);


struct starpu_codelet cl_analyse = {
  .where = STARPU_CPU,
  .cpu_funcs = {_qrm_starpu_analyse_cpu_func, NULL},
  .nbuffers = 1,
  .modes = {STARPU_RW},
  .name = "Analyse"
};

/* submit an solve task
**/
void _qrm_insert_analyse_task_c(qrm_dscr_type *qrm_dscr, void *qrm_mat,
                                starpu_data_handle_t adata_handle,
                                char transp){

  int ret;

  /* printf("> insert analysis, transp:%c  %d\n", transp,qrm_dscr->ctx); */
  ret = starpu_task_insert(&cl_analyse,
                           STARPU_VALUE,     &qrm_dscr, sizeof(void *),
                           STARPU_VALUE,     &qrm_mat, sizeof(void *),
                           STARPU_VALUE,     &transp,  sizeof(char),
                           STARPU_RW,        adata_handle,
                           STARPU_PRIORITY,  1,
                           STARPU_SCHED_CTX, (unsigned)qrm_dscr->ctx,
                           0);
  STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_insert");

  return;

}

/* Many wrappers for the unpack_args routine which is not fortran
   interoperable due the variable length arguments list */

void _qrm_starpu_unpack_args_geqrt(void *cl_arg,
                                   int *bk, int *bi,
                                   int *mbs, int *nbs, int *ib,
                                   int **stair)
{
  starpu_codelet_unpack_args(cl_arg,
                             bk, bi,
                             mbs, nbs, ib,
                             stair);
  return;
}

void _qrm_starpu_unpack_args_gemqrt(void *cl_arg,
                                    int *bk, int *bi, int *bj,
                                    int *mbs, int *nbs, int *ib,
                                    int **stair)
{
  starpu_codelet_unpack_args(cl_arg,
                             bk, bi, bj,
                             mbs, nbs, ib,
                             stair);
  return;
}

void _qrm_starpu_unpack_args_tpqrt(void *cl_arg,
                                   int *bk, int *bi, int *bl,
                                   int *mbs, int *nbs, int *ib, char *ts,
                                   int **stair)
{
  starpu_codelet_unpack_args(cl_arg,
                             bk, bi, bl,
                             mbs, nbs, ib, ts,
                             stair);
  return;
}

void _qrm_starpu_unpack_args_tpmqrt(void *cl_arg,
                                    int *bk, int *bi, int *bl, int *bj,
                                    int *mbs, int *nbs, int *ib, char *ts,
                                    int **stair)
{
  starpu_codelet_unpack_args(cl_arg,
                             bk, bi, bl, bj,
                             mbs, nbs, ib, ts,
                             stair);
  return;
}


/* TODO: qrm_mat should be void **  ??? */
void _qrm_starpu_unpack_args_subtree(void *cl_arg,
                                     void *qrm_dscr,
                                     void *qrm_mat,
                                     int *rootnum)
{
  starpu_codelet_unpack_args(cl_arg,
                             qrm_dscr,
                             qrm_mat,
                             rootnum);
  return;
}

void _qrm_starpu_unpack_args_init(void *cl_arg,
                                  void *qrm_dscr,
                                  void *qrm_mat,
                                  int *fnum)
{
  starpu_codelet_unpack_args(cl_arg,
                             qrm_dscr,
                             qrm_mat,
                             fnum);
  return;
}

void _qrm_starpu_unpack_args_asm(void *cl_arg,
                                 void *qrm_mat, int *fnum,
                                 int *br, int *bc)
{
  starpu_codelet_unpack_args(cl_arg,
                             qrm_mat, fnum,
                             br, bc);
  return;
}


void _qrm_starpu_unpack_args_clean(void *cl_arg,
                                   void *qrm_dscr,
                                   void *qrm_mat,
                                   int *fnum)
{
  starpu_codelet_unpack_args(cl_arg,
                             qrm_dscr,
                             qrm_mat,
                             fnum);
  return;
}

void _qrm_starpu_unpack_args_apply(void *cl_arg,
                                   void *qrm_dscr,
                                   void *qrm_mat,
                                   void **b,
                                   int *inode,
                                   char *transp
                                   )
{
  starpu_codelet_unpack_args(cl_arg,
                             qrm_dscr,
                             qrm_mat,
                             b,
                             inode,
                             transp);
  return;
}

void _qrm_starpu_unpack_args_solve(void *cl_arg,
                                   void *qrm_dscr,
                                   void *qrm_mat,
                                   void **b,
                                   void **x,
                                   int *inode,
                                   char *transp
                                   )
{
  starpu_codelet_unpack_args(cl_arg,
                             qrm_dscr,
                             qrm_mat,
                             b,
                             x,
                             inode,
                             transp);
  return;
}



void _qrm_starpu_unpack_args_analyse(void *cl_arg,
                                     void *qrm_dscr,
                                     void *qrm_mat,
                                     char *transp
                                     )
{
  starpu_codelet_unpack_args(cl_arg,
                             qrm_dscr,
                             qrm_mat,
                             transp);
  return;
}
