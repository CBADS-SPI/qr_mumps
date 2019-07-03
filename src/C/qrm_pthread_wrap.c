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

#include <pthread.h>
#include <stdlib.h>

int pthread_mutex_init_c(void *mutex, void *attr) {

  return pthread_mutex_init((pthread_mutex_t *)mutex, 
                            (pthread_mutexattr_t *)attr);
}

void qrm_alloc_pthread_mutex_c(void **ptr) {
  
  *ptr = (void *)malloc(sizeof(pthread_mutex_t));

  return;
}

int pthread_cond_init_c(void *cond, void *attr) {

  return pthread_cond_init((pthread_cond_t *)cond, 
                            (pthread_condattr_t *)attr);
}

void qrm_alloc_pthread_cond_c(void **ptr) {
  
  *ptr = (void *)malloc(sizeof(pthread_cond_t));

  return;
}

void qrm_dealloc_pthread_cond_c(void **ptr) {
  
  free(*ptr);
  *ptr = NULL;

  return;
}

void qrm_dealloc_pthread_mutex_c(void **ptr) {
  
  free(*ptr);
  *ptr = NULL;

  return;
}


int pthread_cond_wait_c(void *cond,
                        void *mutex) {

  return pthread_cond_wait((pthread_cond_t *)cond,
                           (pthread_mutex_t *) mutex);
}

int pthread_cond_signal_c(void *cond) {
  
  pthread_cond_signal((pthread_cond_t *)cond);
}
