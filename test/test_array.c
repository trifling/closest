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

int random_test( int nd, int ni, int no, int nq ) {
   
   srand(time(NULL));

   double *xi = calloc( nd*ni, sizeof(double) );

   /* generate random points */
   for( int i=0; i<ni; i++ ) {
      for( int k=0; k<nd; k++ ) {
         xi[k+i*nd] = (double)rand()/(double)RAND_MAX; 
      }
   }

   cell_t *cell = cell_init( nd, ni, xi, 0.0 );
   brut_t *brut = brut_init( nd, ni, xi );

   int *i_cell = calloc( nq*no, sizeof(int) );
   int *i_brut = calloc( nq*no, sizeof(int) );
   
   double *d_cell = calloc( nq*no, sizeof(double) );
   double *d_brut = calloc( nq*no, sizeof(double) );

   double *xq = calloc( nd*nq, sizeof(double) );
   for( int i=0; i<nq; i++ ) {
      for( int k=0; k<nd; k++ ) {
         xq[k+i*nd] = (double)rand()/(double)RAND_MAX; 
      }
   }

   printf("Performing %d look-ups for the %d closest point(s), in a set of %d points,"
          " to a random point X in %d-dimensional space...", nq, no, ni, nd );
   fflush(stdout);

   brut_knearest( brut, nq, xq, no, i_brut, d_brut ); 
   cell_knearest( cell, nq, xq, no, i_cell, d_cell ); 

   for( int i=0; i<nq; i++ ) {

      for( int l=0; l<no; l++ ) {
         
         int m = l+i*no;

         if( i_cell[m] != i_brut[m] ) {
            
            fprintf(stderr,"\n\nERROR at iter i=%u!\n", i );
            fprintf(stderr,"X = ( ");
            for( int k=0; k<nd; k++ ) 
               fprintf(stderr," %f ",xq[k]);
            fprintf(stderr,")\n");

            for( int j=0; j<no; j++ ) {
               int m = j+i*no;
               fprintf(stderr, "CELL %02d -> i=%02d d=%f x=(%f,%f)\n", l, i_cell[m], d_cell[m], xi[0+i_cell[m]*nd], xi[1+i_cell[m]*nd] );
               fprintf(stderr, "BRUT %02d -> i=%02d d=%f x=(%f,%f)\n", l, i_brut[m], d_brut[m], xi[0+i_brut[m]*nd], xi[1+i_brut[m]*nd] );
            }

            fprintf(stderr,"Test failed!\n");
            return -1;
         }
      }
   }

   free(xi);
   free(i_cell);
   free(i_brut);
   free(d_cell);
   free(d_brut);
   free(xq);
   cell_free(cell);
   printf("  succeeded!\n");
   return 0;
}

int main( int argc, char *argv[] ) {

   if( argc < 5 ) {
      printf("closest_test nd ni no nq\n");
      printf("  nd - dimensions\n");
      printf("  ni - data points\n");
      printf("  no - number of nearest neighbors to ask for\n");
      printf("  nq - number of queries\n");
      exit(-1);
   }

   int nd = atoi(argv[1]);
   int ni = atoi(argv[2]);
   int no = atoi(argv[3]);
   int nq = atoi(argv[4]);

   return random_test( nd, ni, no, nq );
}

