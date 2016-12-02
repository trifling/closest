# Closest 
 
A C library for the k-nearest neighbor search optimization problem in N-dimensions. 
Three strategies are implemented: brute-force, cull-based and cell-based. 

Closest is small, relatively efficient, and can be easily integrated.


# Code Example

Find the `k`=10 nearest neighbors to a point `x` in an array `data` of 1000 points in `d`=20 dimensions
using the cull, cell and brute-force methods.

```c
   int d = 20;
   int k = 10;
   int n = 1000;

   /* 
      double *data[n*d] has been filled with 1000 points in d=20 dimensions, each
      of the 20 coordinates for a point followed by the next until all 1000 
      are specified.
     
      x is the point where the nearest neighbors query is done
   */

   /* 
    * using the cull method 
    */ 

   /* initialize the search */
   cull_t *cull = cull_init( d, n, data ); 

   int i_cull[k];
   int d_cull[k];
   /* cull_knearest returns 
        in i_cull, the index for the points in the `data` array of the k-nearest neighbors to x */
        in d_cull, the distances for those same points to x */
   cull_knearest( cull, x, k, i_cull, d_cull ); 
   cull_free(cull);

   /* 
    * using the cell method 
    */ 

   /* initialize the search */
   cell_t *cell = cell_init( d, n, -1, data );
   
   int i_cell[k];
   int d_cell[k];
   /* cell_knearest returns 
        in i_cell, the index for the points in the `data` array of the k-nearest neighbors to x */
        in d_cell, the distances for those same points to x */
   cell_knearest( cell, x, k, i_cell, d_cell ); 
   cell_free(cell);

   /* 
    * using the brute-force method 
    */ 
   int i_brut[k];
   int d_brut[k];
   bruteforce_knearest( d, n, data, x, k, i_brut, d_brut ); 

```

# Build and install

Closest uses cmake as buld system. To build and install a dynamic version of closest in /usr/local:

```bash
$ mkdir bld && cd bld && cmake -DCMAKE_INSTALL_PREFIX=/usr/local ../ && make && make install && cd ..
``` 

## GNU GPL v3 License
Copyright (c) 2011-2016, Daniel Pena 

