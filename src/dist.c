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

#include <math.h>

double L1metric( int nd, double *x, double *y, void *data ) {
   (void)data;
   double t = 0.0;
   for( int k=0; k<nd; k++ ) {
      t += fabs(x[k]-y[k]);
   }
   return t;
}

double L2metric( int nd, double *x, double *y, void *data ) {
   (void)data;
   double t = 0.0;
   for( int k=0; k<nd; k++ ) {
      double z = ( x[k] - y[k] );
      t += z*z;
   }
   return sqrt(t);
}

double Lpmetric( int nd, double *x, double *y, void *data ) {
   double f = *(double *)data;
   double t = 0.0;
   for( int k=0; k<nd; k++ ) {
      t += pow( fabs(x[k]-y[k]), f );
   }
   return pow(t,1.0/f);
}

double LMmetric( int nd, double *x, double *y, void *data ) {
   double *M = (double *)data;
   double s = 0.0;
   for( int l=0; l<nd; l++ ) {
      double t = 0.0;
      for( int k=0; k<nd; k++ ) {
         t += M[k+l*nd]*(x[k]-y[k]);
      }
      s += (x[l]-y[l])*t; 
   }
   return s;
}
      
