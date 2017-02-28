/* 
   This file is part of closest, a library for k-nearest neighbors (kNN) search
   in n-dimensions.  
   Copyright (C) 2011-2016 Daniel Pena <trifling.github@gmail.com>

   Closest is free software: you can redistribute it and/or modify it under the
   terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version.

   Closest is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
   FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
   details.

   You should have received a copy of the GNU General Public License along with
   closest. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CLOSEST_AUXF_H
#define CLOSEST_AUXF_H

#include <math.h>
#include <limits.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
  
typedef double (*metric)( int, double *, double *, void * );
double L1metric( int nd, double *x, double *y, void *data ); 
double L2metric( int nd, double *x, double *y, void *data ); 
double Lpmetric( int nd, double *x, double *y, void *data ); 
double LMmetric( int nd, double *x, double *y, void *data ); 

#ifndef MAX
#define MAX(x,y) (((x) > (y) ? (x) : (y)))
#endif
#ifndef MIN
#define MIN(x,y) (((x) < (y) ? (x) : (y)))
#endif

/* RESULT */
typedef struct {
   int i;
   double r;
} result_t;
inline result_t *result_init( int i, double r ) {
   result_t *res = (result_t *)calloc( 1, sizeof(result_t) );
   res->i = i;
   res->r = r;
   return res;
}
inline void result_free( result_t *r ) {
   free(r);
}

/* AUX */
inline double dminval( int n, int stride, double *x ) {
   double r = DBL_MAX;
   int k=0;
   for( int i=0; i<n; i++ ) {
      if( x[k] < r )
         r = x[k];
      k+=stride;
   }
   return(r);
}
inline double dmaxval( int n, int stride, double *x ) {
   double r = -DBL_MAX;
   int k=0;
   for( int i=0; i<n; i++ ) {
      if( x[k] > r )
         r = x[k];
      k+=stride;
   }
   return(r);
}
inline double ddot( int n, double *x, double *y ) {
   double s = 0.0;
   for( int i=0; i<n; i++ ) 
      s += x[i]*y[i];
   return s;
}
inline int udot( int n, int *x, int *y ) {
   int s = 0.0;
   for( int i=0; i<n; i++ ) 
      s += x[i]*y[i];
   return s;
}
inline int any_equal( int n, int *x, int *y ) {
   for( int i=0; i<n; i++ ) 
      if( x[i] == y[i] )
         return 1;
   return 0;
}
inline double L2dist2( int nd, double *x1, double *x2 ) {
   double t = 0.0;
   for( int k=0; k<nd; k++ ) {
      double z = ( x1[k] - x2[k] );
      t += z*z;
   }
   return t;
}
inline void pdvec( const char *s, int n, double *x ) {
   printf("%s = (",s); 
   for( int k=0; k<n; k++ ) {
      printf(" %f", x[k] );
   }
   printf(" ) \n");
}
inline void pivec( const char *s, int n, int *x ) {
   printf("%s = (",s); 
   for( int k=0; k<n; k++ ) {
      printf(" %d", x[k] );
   }
   printf(" ) \n");
}

#endif

