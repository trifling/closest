/* 
   This file is part of closest, a library for cull-based and 
   cell-based spatial search of k-nearest neighbors in n-dimensions.
   Copyright (C) 2011-2016 Daniel Pena <trifling.github@gmail.com>

   Closest is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   Closest is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with closest. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CLOSEST_AUXF_H
#define CLOSEST_AUXF_H

#include <math.h>
#include <limits.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef MAX
#define MAX(x,y) (((x) > (y) ? (x) : (y)))
#endif
#ifndef MIN
#define MIN(x,y) (((x) < (y) ? (x) : (y)))
#endif

/* types */
struct cull_s {
   double *xi;
   int *bmap;
   int *fmap;
   double *oset;
   int ni;
   int nd;
   double *llim;
   double *ulim;
   double *side;
   int *ix;
   int *lst;
   double *lrs;
   char *use;
   double probability;
   double eps;
};

struct cell_s {
   int nd, ni, nr, nc;
   int *stride; 
   int *cell;   
   int *next;   
   double *xi;
   double *xmin;
   double *xmax;
   double *xdel;
};

/* ranking functions */
typedef struct {
   int index;
   double value;
} closest_rank_t;

inline int closest_rank_compare( const void *a, const void *b ) {
   double x = ((closest_rank_t *)a)->value;
   double y = ((closest_rank_t *)b)->value;
   return ((x) < (y) ? -1 : ((x) == (y) ? 0 : 1));
} 

inline void closest_rank( double *data, int nmemb, int *index ) {
   closest_rank_t *elems = calloc( nmemb, sizeof(closest_rank_t) ); 
   for( int i=0; i<nmemb; i++ ) {
      elems[i].index = i;
      elems[i].value = data[index[i]];
   }
   qsort( elems, nmemb, sizeof(elems[0]), closest_rank_compare );
   for( int i=0; i<nmemb; i++ ) {
      index[i] = elems[i].index;
   }
   free(elems);
}

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
inline double L2dist( int nd, double *x1, double *x2 ) {
   return sqrt(L2dist2(nd,x1,x2));
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

