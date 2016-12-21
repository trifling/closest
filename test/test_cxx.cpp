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
#include <closest.h>

int random_test( int nd, int ni, int no, int nr ) {
   
   srand(time(NULL));

   double *xi = new double[nd*ni];

   /* generate random points */
   for( int i=0; i<ni; i++ ) {
      for( int k=0; k<nd; k++ ) {
         xi[k+i*nd] = (double)rand()/(double)RAND_MAX; 
      }
   }

   closest::Cell *cell = new closest::Cell( nd, ni, xi, 0.0 );
   closest::Brut *brut = new closest::Brut( nd, ni, xi );

   int *i_cell = new int[no]; 
   int *i_brut = new int[no];
   
   double *d_cell = new double[no];
   double *d_brut = new double[no];

   double x[nd];

   printf("Performing %d look-ups for the %d closest point(s), in a set of %d points,"
          " to a random point X in %d-dimensional space...", nr, no, ni, nd );
   fflush(stdout);

   for( int i=0; i<nr; i++ ) {

      for( int k=0; k<nd; k++ ) {
         x[k] = (double)rand()/(double)RAND_MAX; 
      }
         
      cell->knearest( 1, x, no, i_cell, d_cell ); 
      brut->knearest( 1, x, no, i_brut, d_brut ); 

      for( int l=0; l<no; l++ ) {
         if( i_cell[l] != i_brut[l] ) {
            
            fprintf(stderr,"\n\nERROR at iter i=%u!\n", i );
            fprintf(stderr,"X = ( ");
            for( int k=0; k<nd; k++ ) 
               fprintf(stderr," %f ",x[k]);
            fprintf(stderr,")\n");

            for( int l=0; l<no; l++ ) {
               fprintf(stderr, "CELL %02d -> i=%02d d=%f x=(%f,%f)\n", l, i_cell[l], d_cell[l], xi[0+i_cell[l]*nd], xi[1+i_cell[l]*nd] );
               fprintf(stderr, "BRUT %02d -> i=%02d d=%f x=(%f,%f)\n", l, i_brut[l], d_brut[l], xi[0+i_brut[l]*nd], xi[1+i_brut[l]*nd] );
            }

            fprintf(stderr,"Test failed!\n");
            return -1;
         }
      }
   }

   delete xi;
   delete i_cell;
   delete i_brut;
   delete d_cell;
   delete d_brut;
   delete cell;
   delete brut;
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

