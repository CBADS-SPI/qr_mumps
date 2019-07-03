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



/*##############################################################################################
** @file secs.c
** FIXME: add comments
**
** @date    05-04-2011
** @author  Alfredo Buttari
** @version 0.0.1
**
**##############################################################################################*/


#include <sys/time.h>

void secs_(double *s)
{
  struct timeval tp;
  struct timezone tzp;
  int i;

  i = gettimeofday(&tp,&tzp);
  *s =  ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );

  return;
}


void usecs_(double *s){
  struct timeval t;
  struct timezone tzp;

  gettimeofday(&t,&tzp);
  *s = (double) t.tv_sec*1000000+t.tv_usec;

  return;

}
