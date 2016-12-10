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

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <unistd.h>
#include <time.h>
#include "auxf.h"

typedef struct {
   int index;
   double value;
} closest_rank_t;


#define RANK_PREFIX double
#define RANK_DATA_TYPE double 

#define HEAP_PREFIX result 
#define HEAP_DATA_TYPE result_t
#define HEAP_CMP(x,y) ((x->r) < (y->r) ? -1 : ((x->r) == (y->r) ? 0 : 1))
#include "maxheap.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define DEFAULT_PROBABILITY 0.3

/* types */
typedef struct {
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
   char *use;
   double prob;
   double eps;
   metric Lfun;
   void *Ldata;
} cull_t;

/* uniform distribution epsilon for a sphere */
/*static double hsphere( int nd, double *ext, int ni, double p ) {*/
   /*double eps = nd*tgamma(nd/2.0)/(2.0*pow(M_PI,nd/2.0))*(1.0-pow(1.0-p,1.0/ni));*/
   /*for( int k=0; k<nd; k++ ) */
      /*eps *= ext[k];*/
   /*return pow(eps,1.0/nd);*/
/*}*/

/* uniform distribution epsilon for a cube */
static double hcube( int nd, double *ext, int ni, double p ) {
   double eps = pow(0.5,nd)*(1.0-pow(1.0-p,1.0/ni));
   for( int k=0; k<nd; k++ ) 
      eps *= ext[k];
   return pow(eps,1.0/nd);
}

void cull_set_metric( cull_t *cull, metric fun, void *data ) {
   cull->Ldata = data;
   if( fun == NULL ) 
      cull->Lfun = L2metric;
   else
      cull->Lfun = fun;
}

/* pre-process the node set to obtain an ordered set, along with a forward and backwards mapping. */
cull_t *cull_init( int nd, int ni, double *xi, double prob ) { 
   
   cull_t *c; 
   c = malloc( sizeof(cull_t) ); 
   c->nd = nd; 
   c->ni = ni; 
   c->xi = xi; 
   c->bmap = calloc( ni, sizeof(int) ); 
   c->fmap = calloc( ni*nd, sizeof(int) );
   c->oset = calloc( ni*nd, sizeof(double) );
   c->llim = calloc( nd, sizeof(double) ); 
   c->ulim = calloc( nd, sizeof(double) ); 
   c->side = calloc( nd, sizeof(double) );

   c->use = calloc( c->ni, sizeof(char) );
   c->ix  = calloc( c->ni, sizeof(int) );
   c->lst = calloc( c->ni, sizeof(int) );

   c->Lfun = L2metric;
   c->Ldata = NULL;

   /* limits */
   for( int k=0; k<nd; k++ ) {
      c->llim[k] = dminval( ni, nd, xi+k );
      c->ulim[k] = dmaxval( ni, nd, xi+k );
      c->side[k] = c->ulim[k] - c->llim[k]; 
   }

   /* create backward and forward maps and ordered sets */
   int *idx = calloc( ni, sizeof(int) );
   double *tmp = calloc( ni, sizeof(double) );
   for( int k=0; k<nd; k++ ) {
      
      /* order taking into account the stride */
      for( int i=0; i<ni; i++ ) { 
         idx[i] = i;
         tmp[i] = c->xi[k+i*nd];
      }

      double_rank( idx, tmp, ni ); 
      /*closest_rank( xi+k, ni, idx ); */

      for( int i=0; i<ni; i++ ) {
         c->oset[k*ni+i] = c->xi[k+idx[i]*nd];
         c->fmap[k*ni+idx[i]] = i;
         
         /*c->oset[k*ni+i] = c->xi[k+idx[i]*nd];*/
         /*c->fmap[k*ni+idx[i]] = i;*/
      }
      if( k==0 ) {
         for( int i=0; i<ni; i++ ) {
            c->bmap[i] = idx[i]; 
         }
      }
   }
   free( idx );

   /* numerical parameters */
   if( prob <= 0.0 ) 
      c->prob = DEFAULT_PROBABILITY;
   else
      c->prob = prob;
   c->eps = hcube( c->nd, c->side, c->ni, c->prob );
   return( c );
}


void cull_free( cull_t *c ) {
   free( c->bmap );
   free( c->fmap );
   free( c->oset );
   free( c->llim );
   free( c->ulim );
   free( c->side );
   free( c->use  );
   free( c->ix   );
   free( c->lst  );
   free( c );
}


/* bin search with lower bias */ 
static inline int bbinsearch( size_t n, double *ordered, double val ) {
   int b = 0;
   int t = n-1;

   if( val > ordered[n-1] )
      return n-1;

   while( t > b+1 ) {
      int c = (b + t)  / 2;
      if( val < ordered[c] ) 
         t = c;
      else
         b = c;
   }
   return b;
}

/* bin search with upper bias */ 
static inline int tbinsearch( size_t n, double *ordered, double val ) {
   int b = 0;
   int t = n-1;

   if( val < ordered[0] )
      return 0;

   while( t > b+1 ) {
      int c = (b + t)  / 2;
      if( val > ordered[c] ) 
         b = c;
      else
         t = c;
   }
   return t;
}

static inline int hyperrect( cull_t *c, double *x, double eps, int *lst ) {

   /* l - index from the ordered set to original set (bmap) */
   /* i - index in the original set */
   /* j - index from original set to ordered set */
   /* k - index in dimensions */
   /* m, n, o - running indexes */
   
   /* initial culling along dim 0, add candidates to the list */
   int bot = bbinsearch( c->ni, c->oset, x[0]-eps );
   int top = tbinsearch( c->ni, c->oset, x[0]+eps );

   int n = 0;
   for( int l=bot; l<=top; l++ ) { 
      int i = c->bmap[l]; 
      if( c->use[i] == 0 )  
         lst[n++] = c->bmap[l];
   }

   /* trim list with binary searches along the other dims */
   for( int k=1; k<c->nd; k++ ) {

      bot = bbinsearch( c->ni, c->oset+k*c->ni, x[k]-eps );
      top = tbinsearch( c->ni, c->oset+k*c->ni, x[k]+eps );

      int m = n; n = 0;
      for( int o=0; o<m; o++ ) {
         int i = lst[o];
         int j = c->fmap[k*c->ni+i];
         if( j < bot || j > top ) 
            continue;

         lst[n++] = i;
      }
   }

   /* mark points already used */
   for( int l=0; l<n; l++ )  
      c->use[lst[l]] = 1;

   return n;
}


int cull_knearest( cull_t *c, int nq, double *xq, int no, int *index, double *distance ) {

   /* epsilon estimation */
   double eps = c->eps * pow((double)no,1.0/(double)c->nd);
   result_heap_t *heap = result_heap_init( no, -1 );
  
   for( int iq=0; iq<nq; iq++ ) {

      result_heap_reset( heap );

      double *x = xq + iq*c->nd;

      /* reset points already used array */
      for( int i=0; i<c->ni; i++ ) 
         c->use[i] = 0;
      
      double rmax = DBL_MAX;
      int n = 0;

      do {
         
         int m = hyperrect( c, x, eps, c->lst+n );

         for( int j=n; j<n+m; j++ ) {
            int i = c->lst[j];
            double r = c->Lfun( c->nd, x, c->xi+i*c->nd, c->Ldata ); 
            if( r < rmax ) {
               result_heap_insert( heap, result_init(i,r) );
               if( heap->len == no )
                  rmax = result_heap_peek( heap )->r; 
            }
         }
         n += m;
      
         if(  n >= no )
            break;
         
         eps *= 1.2;

      } while(1); 

      if( rmax > eps ) 
         eps = rmax;
      else  
         eps = sqrt((double)c->nd)*eps; 
       
      int m = hyperrect( c, x, eps, c->lst+n );
      for( int j=n; j<n+m; j++ ) {
         int i = c->lst[j];
         double r = c->Lfun( c->nd, x, c->xi+i*c->nd, c->Ldata ); 
         if( r < rmax ) {
            result_heap_insert( heap, result_init(i,r) );
            if( heap->len == no )
               rmax = result_heap_peek( heap )->r; 
         }
      }
      n += m;

      /* sort the list according to metric */
      for( int i=no-1; i>=0; i-- ) {
         result_t *e = result_heap_pop( heap );
         if(e) {
            index[i+iq*no] = e->i;
            distance[i+iq*no] = e->r;
            free(e);
         }
      }
   }

   return no*nq;
}

