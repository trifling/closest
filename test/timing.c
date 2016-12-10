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

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <closest.h>
#include "kdtree.h"

int timing_test( int nd, int ni, int no, int nr, double *t ) {

   clock_t start, end;

   double *xi = calloc( nd*ni, sizeof(double) );
   double *x  = calloc( nd*nr, sizeof(double) );
   
   /* generate random points */
   struct drand48_data rd;
   srand(time(NULL));
   short unsigned int seed16v[3];
   seed16v[0] = rand();
   seed16v[1] = rand();
   seed16v[2] = rand(); 
   seed48_r( seed16v, &rd );
   for( int i=0; i<ni*nd; i++ ) 
      drand48_r( &rd, xi+i );
   for( int i=0; i<nr*nd; i++ ) 
      drand48_r( &rd, x +i );

   double *dst = calloc( no, sizeof(double) );
   int *idx = calloc( no, sizeof(int) );
    
   /* time cull method */
   start = clock();
   cull_t *cull = cull_init( nd, ni, xi ); 
   for( int i=0; i<nr; i++ ) {
      cull_knearest( cull, x+i*nd, no, idx, dst ); 
   }
   cull_free(cull);
   end = clock();
   t[0] = ((double) (end - start)) / CLOCKS_PER_SEC;

   /* time cell method */
   start = clock();
   cell_t *cell = cell_init( nd, ni, -1, xi ); 
   for( int i=0; i<nr; i++ ) {
      cell_knearest( cell, x+i*nd, no, idx, dst ); 
   }
   cell_free(cell);
   end = clock();
   t[1] = ((double) (end - start)) / CLOCKS_PER_SEC;

   /* time bruteforce method */
   start = clock();
   for( int i=0; i<nr; i++ ) {
      bruteforce_knearest( nd, ni, xi, x+i*nd, no, idx, dst ); 
   }
   end = clock();
   t[2] = ((double) (end - start)) / CLOCKS_PER_SEC;

   /* the kdtree method */
   start = clock();
   void *kd = kd_create(nd);
   for( int i=0; i<ni; i++ ) 
      kd_insert( kd, xi+i*nd, NULL );
   for( int i=0; i<nr; i++ )  
      kd_nearest( kd, x+i*nd ); 
   end = clock();
   t[3] = ((double) (end - start)) / CLOCKS_PER_SEC;

   free(xi);
   free(dst);
   free(idx);
   free(x);
   return 0;
}

int main( int argc, char *argv[] ) {
   if( argc < 5 ) {
      printf("closest_timing nd ni no nr\n");
      printf("  nd   - dimensions\n");
      printf("  ni   - data points\n");
      printf("  no   - number of nearest neighbors to ask for\n");
      printf("  nr   - number of times\n");
      exit(-1);
   }

   int nd = atoi(argv[1]);
   int ni = atoi(argv[2]);
   int no = atoi(argv[3]);
   int nr = atoi(argv[4]);

   double times[4];
   timing_test( nd, ni, no, nr, times );

   printf("\n");
   printf("CULL = %f s\n", times[0] );
   printf("CELL = %f s\n", times[1] );
   printf("BRUT = %f s\n", times[2] );
   printf("TREE = %f s\n", times[3] );

   printf("\n");
   printf("CULL and CELL compared to bruteforce\n");
   printf("BRUT / CULL = %f\n", times[2]/times[0] );
   printf("BRUT / CELL = %f\n", times[2]/times[1] );
   
   printf("\n");
   printf("CULL, CELL and BRUT compared to kdtree\n");
   printf("TREE / CULL = %f\n", times[3]/times[0] );
   printf("TREE / CELL = %f\n", times[3]/times[1] );
   printf("TREE / BRUT = %f\n", times[3]/times[2] );
   return 0;
}

