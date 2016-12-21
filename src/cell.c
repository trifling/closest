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

/* CELL */
typedef struct {
   int nd, ni, nr, nc;
   int *stride; 
   int *next; /* ni elems */
   int *cell; /* nr**nd elems */
   double *xi;
   double *xmin;
   double *xmax;
   double *xdel;
   metric Lfun;
   void *Ldata;
} cell_t;

#define HEAP_PREFIX result 
#define HEAP_DATA_TYPE result_t
#define HEAP_CMP(x,y) ((x->r) < (y->r) ? -1 : ((x->r) == (y->r) ? 0 : 1))
#include "maxheap.h"

void cell_set_metric( cell_t *cell, metric fun, void *data ) {
   cell->Ldata = data;
   if( fun == NULL ) 
      cell->Lfun = L2metric;
   else
      cell->Lfun = fun;
}

cell_t *cell_init( int nd, int ni, double *xi, int ncpd ) {

   cell_t *c = calloc( 1, sizeof(cell_t) );
   c->nd = nd;
   c->ni = ni;
   c->xi = xi;
   c->Lfun = L2metric;
   c->Ldata = NULL;

   if( ncpd > 0 )  
      c->nr = ncpd; 
   else
      c->nr = pow( (double)ni/3.0, 1.0/(double)nd ); 

   c->stride = calloc( nd, sizeof(int) );
   c->stride[0] = 1;
   for( int k=1; k<nd; k++ )
      c->stride[k] = c->stride[k-1]*c->nr;
   c->nc = (int) pow((double)c->nr,(double)c->nd);
   c->cell = calloc( c->nc, sizeof(int) );

   c->next = calloc( ni, sizeof(int) );
   for( int l=0; l<c->nc; l++ ) 
      c->cell[l] = 0;
   for( int l=0; l<c->ni; l++ )  
      c->next[l] = 0;

   c->xmin = calloc( nd, sizeof(double) );
   c->xmax = calloc( nd, sizeof(double) );
   c->xdel = calloc( nd, sizeof(double) );
   for( int k=0; k<nd; k++ ) {
      c->xmin[k] = dminval( ni, nd, xi+k );
      c->xmax[k] = dmaxval( ni, nd, xi+k );
      c->xdel[k] = (c->xmax[k] - c->xmin[k]) / (double)c->nr;

      /* error-> cell grid dimensions must be positive */
      if( c->xdel[k] == 0.0 ) {
         free(c->xmin);
         free(c->xmax);
         free(c->xdel);
         free(c);
         return NULL;
      }
   }

   /* Loop on nodes, storing indices in cell and next The loop goes from
    * [1..NI] to be able to avoid the 0 index that cannot be negated to
    * indicate that the node has been marked, thus each time this number is
    * used as an index to an array, for example the original array of nodes xi
    * or next it needs to be decremented by one */
   for( int i=1; i<=ni; i++ ) {

      // Find out to which cell the i-th point belongs to
      int icidx[c->nd];
      for( int k=0; k<c->nd; k++ )  
         icidx[k] = MIN( c->nr-1, (int) ((xi[k+(i-1)*nd]-c->xmin[k])/c->xdel[k]) );
       
      /* find out the cell's linear index */
      int icell = udot( c->nd, icidx, c->stride );

      /* store cell's info */
      
      /* get the head of the list */
      int l = c->cell[icell];
      
      /* insert the i-th node at the head of the list */
      c->cell[icell] = i;
      
      /* push back the list */
      
      /* empty cell -> the first element is the i-th node, mark the end of the list setting next to i */
      if( l == 0 )  
         c->next[i-1] = i;
         
      /* there is something already in the cell, now the next node was the previous head of the list */
      else 
         c->next[i-1] = l;
   }
   return c;
}

void cell_free( cell_t *c ) {
   free( c->cell );
   free( c->next );
   free( c->xmin );
   free( c->xmax );
   free( c->xdel );
   free( c->stride );
   free( c );
   c = NULL;
}

int cell_knearest( cell_t *c, int nq, double *xq, int no, int *index, double *distance ) {

   no = MIN( no, c->ni );

   double delx[c->nd];
   int icidx[c->nd];
   int ilayermin[c->nd];
   int ilayermax[c->nd];
   int irangemin[c->nd];
   int irangemax[c->nd];

   result_heap_t *heap = result_heap_init( no, -1 );

   for( int iq=0; iq<nq; iq++ ) {
      
   double *x = xq + iq*c->nd;
   for( int k=0; k<c->nd; k++ ) 
      delx[k] = x[k] - c->xmin[k];

   for( int k=0; k<c->nd; k++ ) {
      icidx[k] = MIN( c->nr-1, (int) (delx[k]/c->xdel[k]) );
      ilayermin[k] = icidx[k]; // cell indices of the layer whose intersection with the range
      ilayermax[k] = icidx[k]; // defined by [irangemin,irangemax] is currently being searched
      irangemin[k] = 0;
      irangemax[k] = c->nr-1;
   }
   
   result_heap_reset( heap );
   for( int k=0; k<c->nc; k++ )
      c->cell[k] = abs(c->cell[k]);
   int restricted = 0;
   double rmax = DBL_MAX;
   
   /* loop over layers */
   while(1) {

      for( int k=0; k<c->nd; k++ )
         icidx[k] = MAX( irangemin[k], ilayermin[k] );

      /* loop over cells */
      while(1) {

         /* check if current cell index icidx within defined outer cell layer ilayermin, ilayermax */
         if( any_equal(c->nd,icidx,ilayermin) || any_equal(c->nd,icidx,ilayermax) ) {
         
            /* get the cell's linear index */
            int icell = udot( c->nd, icidx, c->stride );

            /* search cell icell for unmarked nodes (if cell is not empty) */
            int l = c->cell[icell];
            
            /* at least one node in cell -> get all nodes in cell */
            if( l > 0 ) {
               
               /* loop on nodes in cell icidx */ 
               while(1) {
                  int ln = c->next[l-1];

                  double r = c->Lfun( c->nd, x, c->xi+(l-1)*c->nd, c->Ldata ); 
                  
                  if( r < rmax ) {
                     result_heap_insert( heap, result_init(l-1,r) );
                     if( heap->len == no )
                        rmax = result_heap_peek( heap )->r; 
                  }

                  if( ln == l ) 
                     break;

                  l = ln;
               }

               /* target no achieved, restrict search */
               if( heap->len >= no && !restricted ) {
                  for( int k=0; k<c->nd; k++ ) {
                     irangemin[k] = MAX(0,MIN(c->nr-1, (int)((delx[k]-rmax)/c->xdel[k]) ));
                     irangemax[k] = MAX(0,MIN(c->nr-1, (int)((delx[k]+rmax)/c->xdel[k]) ));
                  }
                  restricted = 1;
               } 

               /* mark cell as visited */
               c->cell[icell] = -c->cell[icell];

            } /* endif cell is !empty */

         } /* endif any layermin layermax */

         /* Nested cell do loop indexes for each dimension from max(irangemin,ilayermin) up to min(irangemax,ilayermax) */
         int cycle_cells = 0;
         for( int k=0; k<c->nd; k++ ) {
            
            /* if still in range for the current K dimension, increase and break the loop to check for nodes */
            if( icidx[k] < MIN( ilayermax[k],irangemax[k] ) ) {
               icidx[k] += 1;
               cycle_cells = 1;
               break;

            /* otherwise reset the current dimension to origin */
            } else {
               icidx[k] = MAX( irangemin[k], ilayermin[k] );
            }
         }

         if( cycle_cells )
            continue;

         /* All the the nested loops are over, break from cells loop*/
         break;

      } /* end loop over cells */

      /* check if there is still the layer can be pushed further out inside the range*/
      int cycle = 0;
      for( int k=0; k<c->nd; k++ ) {
      
         if( ilayermin[k] > irangemin[k] ) {
            cycle = 1;
            ilayermin[k] -= 1;
         }
         if( ilayermax[k] < irangemax[k] ) {
            cycle = 1;
            ilayermax[k] += 1;
         }
      }
      if( cycle )
         continue;
      
      /* There are no more layers within range, search is over, exit loop */
      break;

   } /* end loop over layers */
   for( int i=no-1; i>=0; i-- ) {
      result_t *e = result_heap_pop( heap );
      if(e) {
         index[i+iq*no] = e->i;
         distance[i+iq*no] = e->r;
         free(e);
      }
   }

   } /* end loop over iq */
   result_heap_free(heap);
   return( no*nq );
}

