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

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <closest.h>

#include "kdtree.h"

int timing_test( int nd, int ni, int no, int nq, double *t ) {

   clock_t start, end;

   double *xi = (double *)calloc( nd*ni, sizeof(double) );
   double *x  = (double *)calloc( nd*nq, sizeof(double) );
   
   /* generate random points */
   srand(time(NULL));
   for( int i=0; i<ni*nd; i++ ) 
      xi[i] = (double)rand()/(double)RAND_MAX; 
   for( int i=0; i<nq*nd; i++ ) 
      x[i] = (double)rand()/(double)RAND_MAX; 

   double *dst = (double *)calloc( no*nq, sizeof(double) );
   int *idx = (int *)calloc( no*nq, sizeof(int) );

   /* time cull method */
   start = clock();
   cull_t *cull = cull_init( nd, ni, xi, 0 ); 
   cull_knearest( cull, nq, x, no, idx, dst ); 
   cull_free(cull);
   end = clock();
   t[0] = ((double) (end - start)) / (double) CLOCKS_PER_SEC;

   /* time cell method */
   start = clock();
   cell_t *cell = cell_init( nd, ni, xi, 0 ); 
   cell_knearest( cell, nq, x, no, idx, dst ); 
   cell_free(cell);
   end = clock();
   t[1] = ((double) (end - start)) / (double) CLOCKS_PER_SEC;

   /* time bruteforce method */
   start = clock();
   brut_t *brut = brut_init( nd, ni, xi );
   brut_knearest( brut, nq, x, no, idx, dst ); 
   brut_free(brut);
   end = clock();
   t[2] = ((double) (end - start)) / (double) CLOCKS_PER_SEC;

   /* time tree method */
   start = clock();
   tree_t *tree = tree_init( nd, ni, xi, 0 ); 
   tree_knearest( tree, nq, x, no, idx, dst ); 
   tree_free(tree);
   end = clock();
   t[3] = ((double) (end - start)) / (double) CLOCKS_PER_SEC;

   /* the kdtree method */
   start = clock();
   void *kd = kd_create(nd);
   for( int i=0; i<ni; i++ ) 
      kd_insert( kd, xi+i*nd, NULL );
   for( int i=0; i<nq; i++ )  
      kd_nearest( kd, x+i*nd ); 
   end = clock();
   t[4] = ((double) (end - start)) / (double) CLOCKS_PER_SEC;

   free(xi);
   free(dst);
   free(idx);
   free(x);
   return 0;
}

int main( int argc, char *argv[] ) {
   if( argc < 5 ) {
      printf("closest_timing nd ni no nq\n");
      printf("  nd   - dimensions\n");
      printf("  ni   - data points\n");
      printf("  no   - number of nearest neighbors to ask for\n");
      printf("  nq   - number of times\n");
      exit(-1);
   }

   int nd = atoi(argv[1]);
   int ni = atoi(argv[2]);
   int no = atoi(argv[3]);
   int nq = atoi(argv[4]);

   double times[10];
   timing_test( nd, ni, no, nq, times );

   printf("%7d %7d %7d %7d %15f %15f %15f %15f %15f\n", nd, ni, no, nq, times[0], times[1], times[2], times[3], times[4] );
   return 0;
}

