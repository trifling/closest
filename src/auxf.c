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

#include "auxf.h"

extern double dminval( int n, int stride, double *x );
extern double dmaxval( int n, int stride, double *x );
extern double ddot( int n, double *x, double *y );
extern int udot( int n, int *x, int *y );
extern int any_equal( int n, int *x, int *y );
extern double L2dist2( int nd, double *x1, double *x2 );
extern double L2dist( int nd, double *x1, double *x2 );
extern void pdvec( const char *s, int n, double *x );
extern void pivec( const char *s, int n, int *x );
extern int closest_rank_compare( const void *a, const void *b );
extern void closest_rank( double *data, int nmemb, int *index );
