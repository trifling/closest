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

#include <stdlib.h>

#include "closest.h"
#include "auxf.h"

int bruteforce_knearest( int nd, int ni, double *xi, double *x, int no, int *index, double *distance ) {

   double *dst = calloc( ni, sizeof(double) );

   for( int i=0; i<ni; i++ ) 
      dst[i] = L2dist( nd, x, xi+i*nd );

   int *idx = calloc( ni, sizeof(int) );
   for( int i=0; i<ni; i++ ) 
      idx[i] = i;
   closest_rank( dst, ni, idx );

   for( int i=0; i<no; i++ ) {
      index[i] = idx[i];
      distance[i] = dst[idx[i]];
   }
   free(idx);
   free(dst);
   return ni;
}

