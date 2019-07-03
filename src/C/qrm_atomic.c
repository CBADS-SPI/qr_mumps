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

#include <stdint.h>

void qrm_atomic_max_int64_t(int64_t *ptr, int64_t value) {
   int64_t old, next;

  while (1) {
    old = *ptr;
    next = ((old)<(value)?(value):(old));
    if (__sync_val_compare_and_swap (ptr, old, next) == old) break;
  };
  
  return;
}


void qrm_atomic_add_int64_t(int64_t *ptr, int64_t value) {

  __sync_add_and_fetch(ptr, value);

}
