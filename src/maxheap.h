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

#if !defined(HEAP_PREFIX) && !defined(RANK_PREFIX)
#error "At least one of HEAP_PREFIX or RANK_PREFIX needs to be set for binary max-heap"
#else

#define MKNAMEAUX(TYPE1,EXT)  TYPE1 ## _ ## EXT 
#define MKNAME(TYPE1,EXT)     MKNAMEAUX(TYPE1,EXT) 

/* 
 * HEAP IMPL
 */
#if defined(HEAP_PREFIX) 
#if !defined(HEAP_DATA_TYPE)
#error "HEAP_PREFIX is set, but the data type is not (HEAP_DATA_TYPE is undefined)"
#endif
#ifndef HEAP_CMP
#define HEAP_CMP(x,y)  ((x) < (y) ? -1 : ((x) == (y) ? 0 : 1))
#endif
#define HEAP_SWAP(x,y) {HEAP_DATA_TYPE *__HEAP_SWAP_t = x; x = y; y = __HEAP_SWAP_t;}
#define HEAP_SHIFT_DOWN    MKNAME(HEAP_PREFIX,heap_shift_down)
#define HEAP_SHIFT_UP      MKNAME(HEAP_PREFIX,heap_shift_up)
#define HEAP_HEAPIFY       MKNAME(HEAP_PREFIX,heap_heapify)
#define HEAP_T             MKNAME(HEAP_PREFIX,heap_t)
#define HEAP_INIT          MKNAME(HEAP_PREFIX,heap_init)
#define HEAP_FREE          MKNAME(HEAP_PREFIX,heap_free)
#define HEAP_POP           MKNAME(HEAP_PREFIX,heap_pop)
#define HEAP_INSERT        MKNAME(HEAP_PREFIX,heap_insert)
#define HEAP_RESET         MKNAME(HEAP_PREFIX,heap_reset)
#define HEAP_PEEK          MKNAME(HEAP_PREFIX,heap_peek)

static __inline void HEAP_SHIFT_DOWN( HEAP_DATA_TYPE **dst, const int start, const int end) {
   int root = start;
   while( (root << 1) <= end ) {
      int child = root << 1;

      if ( (child < end) && (HEAP_CMP( dst[child], dst[child+1] )) < 0)
         child++;

      if (HEAP_CMP( dst[root], dst[child] ) >=0 )
        return;

      HEAP_SWAP( dst[root], dst[child] );
      root = child;
   }
}

static __inline void HEAP_SHIFT_UP( HEAP_DATA_TYPE **dst, const int start, const int end) {
   int child = end;
   while( (child >> 1) >= start ) {
      int root = child >> 1;

      if (HEAP_CMP( dst[root], dst[child] ) >=0 )
        return;

      HEAP_SWAP( dst[root], dst[child] );
      child = root;
   }
}

static __inline void HEAP_HEAPIFY( HEAP_DATA_TYPE **dst, const int size) {
   int start = size >> 1;
   while (start >= 0) {
      HEAP_SHIFT_DOWN( dst, start, size - 1);
      start--;
   }
}

typedef struct {
   int *idx; 
   int len; 
   int max_len; 
   int grw_len; 
   HEAP_DATA_TYPE **dat;
} HEAP_T; 

static __inline HEAP_T *HEAP_INIT( int ini_len, int grow_len ) {
   HEAP_T *heap = (HEAP_T *)calloc( 1, sizeof(HEAP_T) );
   heap->max_len = ini_len;
   heap->grw_len = grow_len;
   heap->dat = (HEAP_DATA_TYPE **)calloc( 1+ini_len, sizeof(HEAP_DATA_TYPE *) );
   heap->len = 0;
   return heap;
}

static __inline void HEAP_RESET( HEAP_T *heap ) {
   for( int i=0; i<heap->len; i++ ) {
      free( heap->dat[i] );
   }
   heap->len = 0;
}

static __inline void HEAP_FREE( HEAP_T *heap ) {
   HEAP_RESET( heap );
   free( heap->dat );
   free( heap );
}

static __inline HEAP_DATA_TYPE *HEAP_POP( HEAP_T *heap ) {
   if( heap->len <= 0 ) 
      return(NULL);
   HEAP_SWAP( heap->dat[0], heap->dat[heap->len-1] );
   heap->len--;
   HEAP_SHIFT_DOWN( heap->dat, 0, heap->len-1 );
   return heap->dat[heap->len];
}

static __inline void HEAP_INSERT( HEAP_T *heap, HEAP_DATA_TYPE *val ) {
   if( heap->len == heap->max_len && heap->grw_len <= 0 ) {
      if( !HEAP_CMP(heap->dat[0],val) ) {
         free(val);
         return;
      }
      heap->dat[heap->len] = val;
      heap->len++;
      HEAP_SHIFT_UP( heap->dat, 0, heap->len-1 );
      free( HEAP_POP( heap ) );
      return;

   } else if( heap->len == heap->max_len ) {
      heap->max_len += heap->grw_len;
      heap->dat = (HEAP_DATA_TYPE **)realloc( heap->dat, 1+heap->max_len*sizeof(HEAP_DATA_TYPE *) );
   } 
   heap->dat[heap->len] = val;
   heap->len++;
   HEAP_SHIFT_UP( heap->dat, 0, heap->len-1 );
   return;
}


static __inline HEAP_DATA_TYPE *HEAP_PEEK( HEAP_T *heap ) {
   return heap->dat[0];
}

#undef HEAP_DATA_TYPE
#undef HEAP_CMP
#undef HEAP_SWAP
#undef HEAP_SHIFT_DOWN 
#undef HEAP_SHIFT_UP   
#undef HEAP_HEAPIFY    
#undef HEAP_T          
#undef HEAP_INIT       
#undef HEAP_FREE       
#undef HEAP_POP        
#undef HEAP_INSERT     
#undef HEAP_RESET   
#undef HEAP_PEEK
#undef HEAP_PREFIX
#endif // heap_impl


/* 
 * RANKING IMPL
 */

#if defined(RANK_PREFIX) 
#if !defined(RANK_DATA_TYPE)
#error "RANK_PREFIX is set, but the data type is not (RANK_DATA_TYPE is undefined)"
#endif
#ifndef RANK_CMP
#define RANK_CMP(x,y)  ((x) < (y) ? -1 : ((x) == (y) ? 0 : 1))
#endif
#define RANK_SWAP(x,y) {int __RANK_SWAP_t = (x); (x) = (y); (y) = __RANK_SWAP_t;}
#define RANK_SHIFT_DOWN    MKNAME(RANK_DATA_TYPE,rank_shift_down)
#define RANK_HEAPIFY       MKNAME(RANK_DATA_TYPE,rank_heapify)
#define RANK               MKNAME(RANK_DATA_TYPE,rank)

static __inline void RANK_SHIFT_DOWN( int *idx, RANK_DATA_TYPE *dst, const int start, const int end) {
   int root = start;
   while( (root << 1) <= end ) {
      int child = root << 1;

      if ( (child < end) && (RANK_CMP( dst[idx[child]], dst[idx[child+1]] )) < 0)
         child++;

      if (RANK_CMP( dst[idx[root]], dst[idx[child]] ) >=0 )
        return;

      RANK_SWAP( idx[root], idx[child] );
      root = child;
   }
}

static __inline void RANK_HEAPIFY( int *idx, RANK_DATA_TYPE *dst, const int size) {
   int start = size >> 1;
   while (start >= 0) {
      RANK_SHIFT_DOWN(idx, dst, start, size - 1);
      start--;
   }
}

static __inline void RANK( int *idx, RANK_DATA_TYPE *dst, const int  size) {
   if (size <= 1) {
      idx[0] = 0;
      return;
   }
   int end = size - 1;
   RANK_HEAPIFY(idx, dst, size);

   while (end > 0) {
      RANK_SWAP( idx[0], idx[end] );
      RANK_SHIFT_DOWN(idx, dst, 0, end - 1);
      end--;
   }
}

#undef RANK_DATA_TYPE
#undef RANK_CMP
#undef RANK_SWAP
#undef RANK_SHIFT_DOWN
#undef RANK_HEAPIFY
#undef RANK
#undef RANK_PREFIX
#endif // RANK impl

#undef MKNAMEAUX
#undef MKNAME

#endif
