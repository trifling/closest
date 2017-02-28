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

#ifndef CLOSEST_H
#define CLOSEST_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct cull_s cull_t;
typedef struct cell_s cell_t;
typedef struct brut_s brut_t;
typedef struct tree_s tree_t;

typedef double (*metric)( int, double *, double *, void * );
double L1metric( int nd, double *x, double *y, void *data ); 
double L2metric( int nd, double *x, double *y, void *data ); 
double Lpmetric( int nd, double *x, double *y, void *data ); 
double LMmetric( int nd, double *x, double *y, void *data ); 

cull_t *cull_init( int nd, int ni, double *xi, double prob );
void cull_set_metric( cull_t *cull, metric fun, void *extra );
int cull_knearest( cull_t *cull, int nq, double *xq, int no, int *index, double *distance );
void cull_free( cull_t *c );

cell_t *cell_init( int nd, int ni, double *xi, int ncpd );
void cell_set_metric( cell_t *cell, metric fun, void *extra );
int cell_knearest( cell_t *cell, int nq, double *xq, int no, int *index, double *distance );
void cell_free( cell_t *c );

brut_t *brut_init( int nd, int ni, double *xi );
void brut_set_metric( brut_t *brut, metric fun, void *extra );
int brut_knearest( brut_t *brut, int nq, double *xq, int no, int *index, double *distance );
void brut_free( brut_t *brut );

tree_t *tree_init( int nd, int ni, double *xi, int ppn );
void tree_set_metric( tree_t *brut, metric fun, void *extra );
int tree_knearest( tree_t *t, int nq, double *xq, int nn, int *index, double *distance );
void tree_free( tree_t *t );

#ifdef __cplusplus
}

namespace closest {

class Cull {
public:
   Cull( int nd, int ni, double *xi, double prob );
   ~Cull();
   int knearest( int nq, double *xq, int no, int *idx, double *dst );
   void set_metric( metric fun, void *data );
private:
   cull_t *m_cull;
   int m_nd, m_ni;
   double *m_xi; 
   double m_prob;
};

class Cell {
public:
   Cell( int nd, int ni, double *xi, int ncpd );
   ~Cell();
   int knearest( int nq, double *xq, int no, int *idx, double *dst );
   void set_metric( metric fun, void *data );
private:
   cell_t *m_cell;
   int m_nd, m_ni;
   double *m_xi; 
   int m_ncpd;
};

class Brut {
public:
   Brut( int nd, int ni, double *xi );
   int knearest( int nq, double *xq, int no, int *idx, double *dst );
   void set_metric( metric fun, void *data );
   ~Brut();
private:
   brut_t *m_brut;
   int m_nd, m_ni;
   double *m_xi; 
};

class Tree {
public:
   Tree( int nd, int ni, double *xi, int ppn );
   ~Tree();
   int knearest( int nq, double *xq, int no, int *idx, double *dst );
   void set_metric( metric fun, void *data );
private:
   tree_t *m_tree;
   int m_nd, m_ni;
   double *m_xi; 
   int m_ppn;
};

}

#endif

#endif /* CLOSEST_H */

