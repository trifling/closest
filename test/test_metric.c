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
#include <float.h>
#include <closest.h>

double infnorm( int nd, double *x, double *y, void *data ) {
   /*double f = *(double *)data;*/
   double t = DBL_MIN;
   for( int k=0; k<nd; k++ ) {
      double a = fabs(x[k]-y[k]);
      if( a > t )
         t = a;
   }
   return t;
}

int random_test( int nd, int ni, int no, int nr ) {
   
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

   int *i_cell = calloc( no, sizeof(int) );
   int *i_brut = calloc( no, sizeof(int) );
   
   double *d_cell = calloc( no, sizeof(double) );
   double *d_brut = calloc( no, sizeof(double) );

   double *x = calloc( nd, sizeof(double) );

   printf("Performing %d look-ups for the %d closest point(s), in a set of %d points,"
          " to a random point X in %d-dimensional space...", nr, no, ni, nd );
   fflush(stdout);

   for( int i=0; i<nr; i++ ) {

      for( int k=0; k<nd; k++ ) {
         x[k+i*nd] = (double)rand()/(double)RAND_MAX; 
      }
      
      brut_set_metric( brut, infnorm, NULL );
      brut_knearest( brut, 1, x, no, i_brut, d_brut ); 

      cell_set_metric( cell, infnorm, NULL );
      cell_knearest( cell, 1, x, no, i_cell, d_cell ); 

      for( int l=0; l<no; l++ ) {
         if(    i_cell[l] != i_brut[l] ) {
            
            fprintf(stderr,"\n\nERROR at iter i=%u!\n", i );
            fprintf(stderr,"X = ( ");
            for( int k=0; k<nd; k++ ) 
               fprintf(stderr," %f ",x[k]);
            fprintf(stderr,")\n");

            for( int l=0; l<no; l++ ) {
               fprintf(stderr, "CELL %02d -> i=%02d d=%f\n", l, i_cell[l], d_cell[l] );
               fprintf(stderr, "BRUT %02d -> i=%02d d=%f\n", l, i_brut[l], d_brut[l] );
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
   free(x);
   cell_free(cell);
   brut_free(brut);
   printf("  succeeded!\n");
   return 0;
}

int main( int argc, char *argv[] ) {

   if( argc < 5 ) {
      printf("closest_test nd ni no nr\n");
      printf("  nd - dimensions\n");
      printf("  ni - data points\n");
      printf("  no - number of nearest neighbors to ask for\n");
      printf("  nr - number of times\n");
      exit(-1);
   }

   int nd = atoi(argv[1]);
   int ni = atoi(argv[2]);
   int no = atoi(argv[3]);
   int nr = atoi(argv[4]);

   return random_test( nd, ni, no, nr );
}

