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

#ifndef CLOSEST_H
#define CLOSEST_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct cull_s cull_t;
typedef struct cell_s cell_t;

cull_t *cull_init( int nd, int ni, double *xi );
void cull_set_probability( cull_t *cull, double probability );
int cull_knearest( cull_t *cull, double *x, int no, int *index, double *distance );
void cull_free( cull_t *cull );

cell_t *cell_init( int nd, int ni, int nr, double *xi );
int cell_knearest( cell_t *c, double *x, int no, int *index, double *distance );
void cell_free( cell_t *c );

int bruteforce_knearest( int nd, int ni, double *xi, double *x, int no, int *index, double *distance );

#ifdef __cplusplus
}

class Closest {
public:
   Closest( int nd, int ni, double *xi );
   ~Closest();
   int cell( double *x, int no, int *idx, double *dst );
   int cull( double *x, int no, int *idx, double *dst );
   int bruteforce( double *x, int no, int *idx, double *dst );
private:
   cell_t *m_cell;
   cull_t *m_cull;
   int m_nd, m_ni;
   double *m_xi; 
};

#endif

#endif /* CLOSEST_H */

