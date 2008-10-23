/*
   Copyright 2005, 2006 Damien Stehl�.

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2 of the License, or (at your
   option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
   more details.

   You should have received a copy of the GNU General Public License along
   with this program; see the file COPYING.  If not, write to the Free
   Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
   02111-1307, USA.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <ctype.h>
#include <limits.h>
#include <gmp.h>

#include "d_mat.h"

double ** d_mat_init(int d, int n)
{
   double ** B;

   B = (double **) malloc (d*sizeof(double*) + n*d*sizeof(double));
   B[0] = (double *) B + d;
	for (long i = 1; i < d; i++) B[i] = B[i-1] + n;

	return B;
}

void d_mat_clear(double ** B)
{
   free(B);
}

void d_mat_print(double ** B, int * expo, int d, int n) 
{
   long i, j; 

   printf("[");
   for (i = 0; i < d; i++) 
   {
      printf("[");
      for (j = 0; j < n; j++) 
      { 
         printf("%E", B[i][j]); 
         if (j < n-1) printf(" "); 
      }
      printf("] * 2^%d\n", expo[i]); 
   }  
   printf("]\n"); 
}

double d_vec_scalar_product(double * vec1, double * vec2, int n)
{
  double sum;

  sum = vec1[0] * vec2[0];
  for (long i = 1; i < n; i++)
     sum += vec1[i] * vec2[i];

  return sum;
} 

double d_vec_norm(double * vec, int n)
{
  double sum;

  sum = vec[0] * vec[0];
  for (long i = 1 ; i < n ; i++)
     sum += vec[i] * vec[i];

  return sum;

} 

