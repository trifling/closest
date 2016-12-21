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

#include <float.h>
#include <limits.h>
#include <math.h>
#include "auxf.h"

#define HEAP_PREFIX result
#define HEAP_DATA_TYPE result_t 
#define HEAP_CMP(x,y) ((x->r) < (y->r) ? -1 : ((x->r) == (y->r) ? 0 : 1))
#include "maxheap.h"

/* NODE */
typedef struct node_s {
   double *xmin, *xmax;  /* bbox */
   int l, u;             /* upper and lower index for the node's content */
   struct node_s *lnode; /* left and right nodes */
   struct node_s *rnode;   
   int sdim;    /* split dimension */
   double sval; /* split value */
   int sidx;    /* split index */

   double rmax;
   double *xc;
} node_t;

/* TREE */
typedef struct {
   /* tree vars */
   int nd, ni;
   double *xi;
   int *idx;
   int nppn; /*n points per leaf node */
   double *xmin, *xmax;
   node_t *root;

   /* search vars */
   int no;
   double *xq;
   double rmax;
   result_heap_t *res;
   metric Lfun;
   void *Ldata;
} tree_t;

/* 
 * AUX FUNCS: ops in indexed arrays 
 */
static inline double didxminval(  const int lower, const int upper, const int *index, const int offset, const int stride, const double *array ) {
   double min = DBL_MAX;
   for( int i=lower; i<= upper; i++ ) {
      int l = index[i];
      double x = array[ offset+l*stride ];
      if( x < min )
         min = x;
   }
   return( min );
}
static inline double didxmaxval(  const int lower, const int upper, const int *index, const int offset, const int stride, const double *array ) {
   double max = DBL_MIN;
   for( int i=lower; i<= upper; i++ ) {
      int l = index[i];
      double x = array[ offset+l*stride ];
      if( x > max )
         max = x;
   }
   return( max );
}
static inline void bbox(  const int l, const int u, const int *idx, const int nd, const double *xi, double *xmin, double *xmax ) {
   for( int k=0; k<nd; k++ )  {
      xmin[k] = didxminval( l, u, idx, k, nd, xi );
      xmax[k] = didxmaxval( l, u, idx, k, nd, xi );
   }
}
static inline void update_bbox( int l, int u, int *idx, int nd, double *xi, int sd, double *ixmin, double *ixmax, double *oxmin, double *oxmax ) {
   for( int k=0; k<nd; k++ ) {
      if( k == sd ) {
         oxmin[k] = didxminval( l, u, idx, k, nd, xi );
         oxmax[k] = didxmaxval( l, u, idx, k, nd, xi );
      } else {
         oxmin[k] = ixmin[k];
         oxmax[k] = ixmax[k];
      }
   }
}

/*
 * NODE
 */
static inline node_t *node_alloc( tree_t *t, int l, int u ) {
   node_t *node = (node_t *)calloc( 1, sizeof(node_t) );
   node->l = l;
   node->u = u;
   node->lnode = NULL;
   node->rnode = NULL; 
   node->xmin = (double *) calloc(t->nd,sizeof(double));
   node->xmax = (double *) calloc(t->nd,sizeof(double));
   return node;
} 

static inline void node_split( const tree_t *t, node_t *node, double *xmin, double *xmax ) {

   /* get the dimension with largest delta into split */
   double maxdel = DBL_MIN;
   int splitd = 0;
   for( int k=0; k<t->nd; k++ ) {
      double d = xmax[k] - xmin[k];
      if( d > maxdel ) {
         splitd = k;
         maxdel = d;
      }
   }
   node->sdim = splitd;
   
   /* splitting method 1 */
   /*node->sval =( xmax[splitd] + xmin[splitd] )*0.5;*/

   /* splitting method 2 */
   /*int i = (l+u)/2;*/
   /*node->sval = t->xi[ splitd + t->idx[i]*t->nd ];*/

   /* average idea taken from kdtree2 */
   int l = node->l;
   int u = node->u;
   double s = 0.0; 
   for( int i=l; i<=u; i++ )  
      s += t->xi[ splitd + t->idx[i]*t->nd ];
   node->sval = s / (double)(u-l+1);

   /* move the node indices along the split dimension left and right the split value */ 
   while( l <= u ) {
      if( t->xi[splitd+t->nd*t->idx[l]] <= node->sval ) {
         l++; 
      } else {
         int tmp = t->idx[l];
         t->idx[l] = t->idx[u];
         t->idx[u] = tmp;
         u--;
      }
   }
   node->sidx = l-1;
} 

node_t *node_init( tree_t *t, int l, int u, double *xmin, double *xmax ) {
   
   node_t *node = node_alloc(t,l,u);
   
   /* if it is a leaf node, calculate the definitive bbox and return */
   if( u-l <= t->nppn ) {
      bbox( l, u, t->idx, t->nd, t->xi, node->xmin, node->xmax );
      
   /* otherwise split it in half and continue */
   } else {
      
      /* split the node in place */
      node_split( t, node, xmin, xmax );

      /* calculate approx bbox for left side, store it temporarily as the node bbox and pass it to the left node init sub */
      update_bbox( l, node->sidx, t->idx, t->nd, t->xi, node->sdim, xmin, xmax, node->xmin, node->xmax );
      node->lnode = node_init( t, l, node->sidx, node->xmin, node->xmax ); 

      /* calculate approx bbox for right side, store it temporarily as the node bbox and pass it to the right node init sub */
      update_bbox( node->sidx+1, u, t->idx, t->nd, t->xi, node->sdim, xmin, xmax, node->xmin, node->xmax );
      node->rnode = node_init( t, node->sidx+1, u, node->xmin, node->xmax );

      /* use children bbox to calculate the node's definitive bbox */
      node->sval = (node->lnode->xmax[node->sdim] + node->rnode->xmin[node->sdim] ) / 2.0; 
      for (int i=0; i<t->nd; i++) {
         node->xmax[i] = MAX(node->lnode->xmax[i],node->rnode->xmax[i]);
         node->xmin[i] = MIN(node->lnode->xmin[i],node->rnode->xmin[i]);
      }
   }

   /* calculate node's center and rmax */
   node->xc = calloc( t->nd, sizeof(double) );
   for( int k=0; k<t->nd; k++ ) {
      node->xc[k] = (node->xmax[k] + node->xmin[k])*0.5;
   }
   node->rmax = t->Lfun( t->nd, node->xc, node->xmax, t->Ldata );

   return(node);
}

static inline void node_free( node_t *node ) {
   if( node->rnode != NULL )
      node_free(node->rnode);
   if( node->lnode != NULL )
      node_free(node->lnode);
   free( node->xmin );
   free( node->xmax );
   free( node );
}

static inline int node_in_range( tree_t *t, node_t *node ) {
   double rt;
   rt = t->Lfun( t->nd, t->xq, node->xc, t->Ldata );
   if( rt < t->rmax + node->rmax )
      return(1);
   return(0);
}

static void node_search( tree_t *t, node_t *node ) { 

   /* if it's an internal node */
   if( node->lnode ) {
      
      /* query val along split dim */
      double x = t->xq[node->sdim]; 

      if( x < node->sval ) {
         node_search( t, node->lnode );
         if( node_in_range( t, node->rnode ) ) 
            node_search( t, node->rnode );
      } else {
         node_search( t, node->rnode );
         if( node_in_range( t, node->lnode ) ) 
            node_search( t, node->lnode );
      }
      return;
   }

   /* otherwise, handle a leaf node */
   for( int i=node->l; i<=node->u; i++ ) {
      
      int l = t->idx[i];
      
      double r = t->Lfun( t->nd, t->xq, t->xi+l*t->nd, t->Ldata );

      /* check if r less than the current max distance */
      if( r > t->rmax ) 
         continue; 
    
      result_heap_insert( t->res, result_init(l,r) );
      if( t->res->len == t->no )  
         t->rmax = result_heap_peek( t->res )->r; 
   } 
}

/*
 * TREE, PUBLIC IFACE
 */

void tree_set_metric( tree_t *tree, metric fun, void *data ) {
   tree->Ldata = data;
   if( fun == NULL ) 
      tree->Lfun = L2metric;
   else
      tree->Lfun = fun;
}

void tree_free( tree_t *t ) {
   node_free( t->root );
   free( t->idx );
   free( t );
}

tree_t *tree_init( int nd, int ni, double *restrict xi, int nppn ) {
   tree_t *t = (tree_t *)calloc( 1, sizeof(tree_t) );
   t->nd = nd;
   t->ni = ni;
   t->xi = xi;
   if( nppn > 0 )
      t->nppn = nppn;
   else
      t->nppn = 32;
   t->root = NULL;
   t->idx = (int *)calloc( ni, sizeof(int) );
   for (int i=0; i<ni; i++) 
      t->idx[i] = i; 
   
   /* data bbox */
   t->xmin = (double *) calloc(t->nd,sizeof(double));
   t->xmax = (double *) calloc(t->nd,sizeof(double));
   bbox( 0, ni-1, t->idx, t->nd, t->xi, t->xmin, t->xmax );

   /* default search params */
   t->Lfun = L2metric;
   t->Ldata = NULL;

   /* build the tree */
   t->root = node_init( t, 0, t->ni-1, t->xmin, t->xmax ); 
   return t;
}


int tree_knearest( tree_t *t, int nq, double *xq, int no, int *index, double *distance ) {
   t->no  = no;
   t->res = result_heap_init( no, -1 );
   for( int iq=0; iq<nq; iq++ ) {
      t->rmax = DBL_MAX;
      t->xq   = xq + iq*t->nd;
   
      node_search( t, t->root );
   
      for( int i=no-1; i>=0; i-- ) {
         result_t *e = result_heap_pop( t->res );
         if(e) {
            index[i+iq*no] = e->i;
            distance[i+iq*no] = e->r;
            free(e);
         }
      }
      result_heap_reset( t->res );
   }
   result_heap_free( t->res );
   return( no*nq );
}

