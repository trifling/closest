
#include "closest.h"

Closest::Closest( int nd, int ni, double  *xi ) 
   : m_cell(nullptr), m_cull(nullptr), m_nd(nd), m_ni(ni), m_xi(xi)
{}       

int Closest::cell( double *x, int no, int *idx, double *dst ) {
   if( m_cell == nullptr ) {
      m_cell = cell_init( m_nd, m_ni, -1, m_xi );
   }
   return cell_knearest( m_cell, x, no, idx, dst );
}

int Closest::cull( double *x, int no, int *idx, double *dst ) {
   if( m_cull == nullptr ) {
      m_cull = cull_init( m_nd, m_ni, m_xi );
   }
   return cull_knearest( m_cull, x, no, idx, dst );
}

int Closest::bruteforce( double *x, int no, int *idx, double *dst ) {
   return bruteforce_knearest( m_nd, m_ni, m_xi, x, no, idx, dst );
}

Closest::~Closest() {
   if( m_cell != nullptr )
      cell_free( m_cell );

   if( m_cull != nullptr )
      cull_free( m_cull );
} 

