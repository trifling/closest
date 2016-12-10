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

#include "auxf.h"

/* BRUT */
typedef struct {
   int nd, ni;
   double *xi;
   metric Lfun;
   void *Ldata;
} brut_t;

#define HEAP_PREFIX result
#define HEAP_DATA_TYPE result_t
#define HEAP_CMP(x,y) ((x->r) < (y->r) ? -1 : ((x->r) == (y->r) ? 0 : 1))
#include "maxheap.h"

void brut_set_metric( brut_t *brut, metric fun, void *data ) {
   brut->Ldata = data;
   if( fun == NULL ) 
      brut->Lfun = L2metric;
   else
      brut->Lfun = fun;
}

brut_t *brut_init( int nd, int ni, double *xi ) {
   brut_t *b = calloc( 1, sizeof(brut_t) );
   b->nd = nd;
   b->ni = ni;
   b->xi = xi;
   b->Lfun = L2metric;
   b->Ldata = NULL;
   return b;
}

int brut_knearest( brut_t *b, int nq, double *xq, int no, int *index, double *distance ) {

   no = MIN( no, b->ni );

   result_heap_t *heap = result_heap_init( no, -1 ); 

   for( int iq=0; iq<nq; iq++ ) {
      double *x = xq + iq*b->nd;

      result_heap_reset( heap );
      double rmax = DBL_MAX;
      for( int i=0; i<b->ni; i++ ) {
         double r = b->Lfun( b->nd, x, b->xi+i*b->nd, b->Ldata );

         if( r > rmax )
            continue;
      
         result_heap_insert( heap, result_init(i,r) );
         if( heap->len == no ) 
            rmax = result_heap_peek( heap )->r; 
      }

      for( int i=no-1; i>=0; i-- ) {
         result_t *e = result_heap_pop( heap );
         if(e) {
            index[i+iq*no] = e->i;
            distance[i+iq*no] = e->r;
            free(e);
         }
      }
   }
   result_heap_free( heap );
   return no*nq;
}

void brut_free( brut_t *b ) {
   free(b);
}

