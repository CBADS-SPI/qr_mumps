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
#ifndef __QRM_STARPU_PROF_H__
#define __QRM_STARPU_PROF_H__


#include <starpu.h>
#include <starpu_profiling.h>
#include <pthread.h>

struct qrm_starpu_prof {

  pthread_mutex_t mutex;
  double activate_time;
  double panel_time;
  double *update_time;
  double assemble_time;
  double clean_time;
  double callback_time;
  double *idle_time;
};

#endif /*__QRM_STARPU_PROF_H__*/
