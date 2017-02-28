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

#include "closest.h"

using namespace closest;

Cull::Cull( int nd, int ni, double *xi, double prob ) 
   : m_cull(nullptr), m_nd(nd), m_ni(ni), m_xi(xi), m_prob(prob) {
   m_cull = cull_init( m_nd, m_ni, m_xi, m_prob );
}       

int Cull::knearest( int nq, double *xq, int no, int *idx, double *dst ) {
   return cull_knearest( m_cull, nq, xq, no, idx, dst );
}

void Cull::set_metric( metric fun, void *data ) {
   cull_set_metric( m_cull, fun, data );
}

Cull::~Cull() {
   cull_free( m_cull );
}

Cell::Cell( int nd, int ni, double *xi, int ncpd ) 
   : m_cell(nullptr), m_nd(nd), m_ni(ni), m_xi(xi), m_ncpd(ncpd) {
   m_cell = cell_init( m_nd, m_ni, m_xi, m_ncpd );
}       

int Cell::knearest( int nq, double *xq, int no, int *idx, double *dst ) {
   return cell_knearest( m_cell, nq, xq, no, idx, dst );
}

void Cell::set_metric( metric fun, void *data ) {
   cell_set_metric( m_cell, fun, data );
}

Cell::~Cell() {
   cell_free( m_cell );
}

Brut::Brut( int nd, int ni, double *xi ) 
   : m_brut(nullptr), m_nd(nd), m_ni(ni), m_xi(xi) {
   m_brut = brut_init( m_nd, m_ni, m_xi );
}       

int Brut::knearest( int nq, double *xq, int no, int *idx, double *dst ) {
   return brut_knearest( m_brut, nq, xq, no, idx, dst );
}

void Brut::set_metric( metric fun, void *data ) {
   brut_set_metric( m_brut, fun, data );
}

Brut::~Brut() {
   brut_free( m_brut );
}

Tree::Tree( int nd, int ni, double *xi, int ppn ) 
   : m_tree(nullptr), m_nd(nd), m_ni(ni), m_xi(xi), m_ppn(ppn) {
   m_tree = tree_init( m_nd, m_ni, m_xi, m_ppn );
}       

int Tree::knearest( int nq, double *xq, int no, int *idx, double *dst ) {
   return tree_knearest( m_tree, nq, xq, no, idx, dst );
}

void Tree::set_metric( metric fun, void *data ) {
   tree_set_metric( m_tree, fun, data );
}

Tree::~Tree() {
   tree_free( m_tree );
}

