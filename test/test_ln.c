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

int random_test( int nd, int ni, int no, int nr, double p ) {
   
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

   int *i_cell_Lf = calloc( no, sizeof(int) );
   int *i_brut_Lf = calloc( no, sizeof(int) );
   int *i_cell_L2 = calloc( no, sizeof(int) );
   int *i_brut_L2 = calloc( no, sizeof(int) );
   
   double *d_cell_Lf = calloc( no, sizeof(double) );
   double *d_brut_Lf = calloc( no, sizeof(double) );
   double *d_cell_L2 = calloc( no, sizeof(double) );
   double *d_brut_L2 = calloc( no, sizeof(double) );

   double *x = calloc( nd, sizeof(double) );

   printf("Performing %d look-ups for the %d closest point(s), in a set of %d points,"
          " to a random point X in %d-dimensional space...", nr, no, ni, nd );
   fflush(stdout);

   for( int i=0; i<nr; i++ ) {

      for( int k=0; k<nd; k++ ) {
         x[k+i*nd] = (double)rand()/(double)RAND_MAX; 
      }
      
      brut_set_metric( brut, Lpmetric, &p );
      brut_knearest( brut, 1, x, no, i_brut_Lf, d_brut_Lf ); 
      brut_set_metric( brut, NULL, NULL );
      brut_knearest( brut, 1, x, no, i_brut_L2, d_brut_L2 ); 

      cell_set_metric( cell, Lpmetric, &p );
      cell_knearest( cell, 1, x, no, i_cell_Lf, d_cell_Lf ); 
      cell_set_metric( cell, NULL, NULL );
      cell_knearest( cell, 1, x, no, i_cell_L2, d_cell_L2 ); 

      for( int l=0; l<no; l++ ) {
         if(    i_cell_Lf[l] != i_brut_Lf[l] 
             || i_cell_L2[l] != i_brut_L2[l] ) { 
            
            fprintf(stderr,"\n\nERROR at iter i=%u!\n", i );
            fprintf(stderr,"X = ( ");
            for( int k=0; k<nd; k++ ) 
               fprintf(stderr," %f ",x[k]);
            fprintf(stderr,")\n");

            for( int l=0; l<no; l++ ) {
               fprintf(stderr, "CELL %02d -> Lf_i=%02d Lf_d=%f L2_i=%02d L2_d=%f \n", l, i_cell_Lf[l], d_cell_Lf[l], i_cell_L2[l], d_cell_L2[l] );
               fprintf(stderr, "BRUT %02d -> Lf_i=%02d Lf_d=%f L2_i=%02d L2_d=%f \n", l, i_brut_Lf[l], d_brut_Lf[l], i_brut_L2[l], d_brut_L2[l] );
            }

            fprintf(stderr,"Test failed!\n");
            return -1;
         }
      }
   }

   free(xi);
   free(i_cell_Lf);
   free(i_brut_Lf);
   free(d_cell_Lf);
   free(d_brut_Lf);
   free(i_cell_L2);
   free(i_brut_L2);
   free(d_cell_L2);
   free(d_brut_L2);
   free(x);
   cell_free(cell);
   printf("  succeeded!\n");
   return 0;
}

int main( int argc, char *argv[] ) {

   if( argc < 6 ) {
      printf("closest_test nd ni no nr\n");
      printf("  nd - dimensions\n");
      printf("  ni - data points\n");
      printf("  no - number of nearest neighbors to ask for\n");
      printf("  nr - number of times\n");
      printf("  p  - Distance power Lp \n");
      exit(-1);
   }

   int nd = atoi(argv[1]);
   int ni = atoi(argv[2]);
   int no = atoi(argv[3]);
   int nr = atoi(argv[4]);
   double p = atof(argv[5]);

   return random_test( nd, ni, no, nr, p );
}

