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

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <unistd.h>
#include <time.h>

#include "closest.h"
#include "auxf.h"

#define INITIAL_PROBABILITY 0.3

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* uniform distribution epsilon for a sphere: for np points it needs multiplying by pow(np,1/nd) */
static double hsphere( int nd, double *ext, int ni, double p ) {
   double eps = nd*tgamma(nd/2.0)/(2.0*pow(M_PI,nd/2.0))*(1.0-pow(1.0-p,1.0/ni));
   for( int k=0; k<nd; k++ ) 
      eps *= ext[k];
   return pow(eps,1.0/nd);
}

/* uniform distribution epsilon for a cube: for np points it needs multiplying by pow(np,1/nd) */
static double hcube( int nd, double *ext, int ni, double p ) {
   double eps = pow(0.5,nd)*(1.0-pow(1.0-p,1.0/ni));
   for( int k=0; k<nd; k++ ) 
      eps *= ext[k];
   return pow(eps,1.0/nd);
}

/* pre-process the node set to obtain an ordered set, along with a forward and backwards mapping. */
cull_t *cull_init(  int nd, int ni, double *xi ) { 
   
   time_t start, end;
   cull_t *cls; 
   cls = malloc( sizeof(cull_t) ); 
   cls->nd = nd; 
   cls->ni = ni; 
   cls->xi = xi; 
   cls->bmap = calloc( ni, sizeof(int) ); 
   cls->fmap = calloc( ni*nd, sizeof(int) );
   cls->oset = calloc( ni*nd, sizeof(double) );
   cls->llim = calloc( nd, sizeof(double) ); 
   cls->ulim = calloc( nd, sizeof(double) ); 
   cls->side = calloc( nd, sizeof(double) );

   cls->lrs = calloc( cls->ni, sizeof(double) );
   cls->use = calloc( cls->ni, sizeof(char) );
   cls->ix  = calloc( cls->ni, sizeof(int) );
   cls->lst = calloc( cls->ni, sizeof(int) );

   /* limits */
   for( int k=0; k<nd; k++ ) {
      cls->llim[k] = dminval( ni, nd, xi+k );
      cls->ulim[k] = dmaxval( ni, nd, xi+k );
      cls->side[k] = cls->ulim[k] - cls->llim[k]; 
   }

   /* create backward and forward maps and ordered sets */
   int *idx = calloc( ni, sizeof(int) );
   for( int k=0; k<nd; k++ ) {
      
      /* order taking into account the stride */
      for( int i=0; i<ni; i++ ) 
         idx[i] = i*nd;
      closest_rank( xi+k, ni, idx ); 

      for( int i=0; i<ni; i++ ) {
         cls->oset[k*ni+i] = cls->xi[k+idx[i]*nd];
         cls->fmap[k*ni+idx[i]] = i;
      }
      if( k==0 ) {
         for( int i=0; i<ni; i++ ) {
            cls->bmap[i] = idx[i]; 
         }
      }
   }
   free( idx );

   /* numerical parameters */
   cull_set_probability( cls, INITIAL_PROBABILITY );
   return( cls );
}

void cull_set_probability( cull_t *cls, double probability ) {
   cls->probability = probability;
   cls->eps = hcube( cls->nd, cls->side, cls->ni, cls->probability );
}

void cull_free( cull_t *cls ) {
   free( cls->bmap );
   free( cls->fmap );
   free( cls->oset );
   free( cls->llim );
   free( cls->ulim );
   free( cls->side );
   free( cls->lrs  );
   free( cls->use  );
   free( cls->ix   );
   free( cls->lst  );
   free( cls );
}


/* bin search with lower bias */ 
static int bbinsearch( size_t n, double *ordered, double val ) {
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
static int tbinsearch( size_t n, double *ordered, double val ) {
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

static int hyperrect( cull_t *cls, double *x, double eps, int *lst ) {

   /* l - index from the ordered set to original set (bmap) */
   /* i - index in the original set */
   /* j - index from original set to ordered set */
   /* k - index in dimensions */
   /* m, n, o - running indexes */
   
   /* initial culling along dim 0, add candidates to the list */
   int bot = bbinsearch( cls->ni, cls->oset, x[0]-eps );
   int top = tbinsearch( cls->ni, cls->oset, x[0]+eps );

   int n = 0;
   for( int l=bot; l<=top; l++ ) { 
      int i = cls->bmap[l]; 
      if( cls->use[i] == 0 )  
         lst[n++] = cls->bmap[l];
   }

   /* trim list with binary searches along the other dims */
   for( int k=1; k<cls->nd; k++ ) {

      bot = bbinsearch( cls->ni, cls->oset+k*cls->ni, x[k]-eps );
      top = tbinsearch( cls->ni, cls->oset+k*cls->ni, x[k]+eps );

      int m = n; n = 0;
      for( int o=0; o<m; o++ ) {
         int i = lst[o];
         int j = cls->fmap[k*cls->ni+i];
         if( j < bot || j > top ) 
            continue;

         lst[n++] = i;
      }
   }

   /* mark points already used */
   for( int l=0; l<n; l++ )  
      cls->use[lst[l]] = 1;
       
   return n;
}


int cull_knearest( cull_t *cls, double *x, int no, int *index, double *distance ) {

   /* epsilon estimation */
   double eps = cls->eps * pow((double)no,1.0/(double)cls->nd);
  
   /* reset points already used array */
   for( int i=0; i<cls->ni; i++ ) 
      cls->use[i] = 0;
   

   double mrs = 0.0;
   int n = 0;

   while( n < no ) {
      int m = hyperrect( cls, x, eps, cls->lst+n );

      eps *= 1.2;

      for( int l=n; l<n+m; l++ ) {
         int i = cls->lst[l];
         cls->lrs[l] = L2dist( cls->nd, x, cls->xi+i*cls->nd ); 
         if( mrs < cls->lrs[l] )
            mrs = cls->lrs[l];
      }
      n += m;
   }

   if( mrs > eps )  
      eps = mrs;
   else  
      eps = sqrt((double)cls->nd)*eps; 
    
   int m = hyperrect( cls, x, eps, cls->lst+n );
   for( int l=n; l<n+m; l++ ) {
      int i = cls->lst[l];
      cls->lrs[l] = L2dist( cls->nd, x, cls->xi+i*cls->nd ); 
   }
   n += m;

   /* sort the list according to metric */
   for( int i=0; i<n; i++ ) 
      cls->ix[i] = i;

   closest_rank( cls->lrs, n, cls->ix ); 

   for( int i=0; i < (n<no?n:no); i++ ) {
      index[i] = cls->lst[cls->ix[i]]; 
      distance[i] = cls->lrs[cls->ix[i]];
   }

   return n;
}

