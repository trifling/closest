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
#include <closest.h>

int random_test( int nd, int ni, int no, int nr ) {

   double *xi = new double[nd*ni];

   /* generate random points */
   struct drand48_data rd;
   srand(time(NULL));
   short unsigned int seed16v[3];
   seed16v[0] = rand();
   seed16v[1] = rand();
   seed16v[2] = rand(); 
   seed48_r( seed16v, &rd );
   for( int i=0; i<ni; i++ ) {
      for( int k=0; k<nd; k++ ) {
         drand48_r( &rd, xi+k+i*nd );
      }
   }

   Closest *cls = new Closest( nd, ni, xi );

   int *i_cull = new int[no];
   int *i_cell = new int[no]; 
   int *i_brut = new int[no];
   
   double *d_cull = new double[no];
   double *d_cell = new double[no];
   double *d_brut = new double[no];

   double x[nd];

   printf("Performing %d look-ups for the %d closest point(s), in a set of %d points,"
          " to a random point X in %d-dimensional space...", nr, no, ni, nd );
   fflush(stdout);

   for( int i=0; i<nr; i++ ) {

      for( int k=0; k<nd; k++ ) {
         double tmp;
         drand48_r( &rd, &tmp );
         x[k] = tmp;
      }
         
      int n_cell = cls->cell( x, no, i_cell, d_cell ); 
      int n_cull = cls->cull( x, no, i_cull, d_cull ); 
      int n_brut = cls->bruteforce( x, no, i_brut, d_brut ); 

      for( int l=0; l<no; l++ ) {
         if( i_cull[l] != i_brut[l] || i_cell[l] != i_brut[l] ) {
            
            fprintf(stderr,"\n\nERROR at iter i=%u!\n", i );
            fprintf(stderr,"X = ( ");
            for( int k=0; k<nd; k++ ) 
               fprintf(stderr," %f ",x[k]);
            fprintf(stderr,")\n");

            for( int l=0; l<no; l++ ) {
               fprintf(stderr, "CELL %02d -> i=%02d d=%f x=(%f,%f)\n", l, i_cell[l], d_cell[l], xi[0+i_cell[l]*nd], xi[1+i_cell[l]*nd] );
               fprintf(stderr, "CULL %02d -> i=%02d d=%f x=(%f,%f)\n", l, i_cull[l], d_cull[l], xi[0+i_cull[l]*nd], xi[1+i_cull[l]*nd] );
               fprintf(stderr, "BRUT %02d -> i=%02d d=%f x=(%f,%f)\n", l, i_brut[l], d_brut[l], xi[0+i_brut[l]*nd], xi[1+i_brut[l]*nd] );
            }

            fprintf(stderr,"Test failed!\n");
            return -1;
         }
      }
   }

   delete xi;
   delete i_cull;
   delete i_cell;
   delete i_brut;
   delete d_cull;
   delete d_cell;
   delete d_brut;
   delete cls;
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
