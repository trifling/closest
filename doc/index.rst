
.. _`GPL-3.0`: https://opensource.org/licenses/GPL-3.0
.. _closest-2.0.zip: https://github.com/trifling/closest/archive/v1.0.zip
.. _github: https://github.com/trifling/closest

This is the documentation for Closest 2.0, last updated December 14, 2016.

Introduction
============
Closest is a C library for k-nearest neighbor search in N-dimensions with 
arbitrary metric functions. 

Closest:
   - is simple to use
   - is flexible (custom distance functions)
   - provides API interfaces for C, C++ and Python
   - written in ISO C99
   - has no external dependencies
   - is small and fast
   - is well documented 

Closest is licensed under the `GPL-3.0`_ (see the LICENSE file for details).

Description
===========

Closest implements four different strategies for the k-nearest neighbor search:
a cell-based method, a kd-tree method, a dimensional culling method and a
relatively efficient brute-force method.

The cell, kd-tree and culling methods work on the assumption that the
underlying set is a vector space with a metric function (not necessarily a
norm). The brute-force method only assumption is that there is a distance
function that is defined in the working set, and in that regard is more general
that the other three.

The default metric is the Euclidean norm (or L2), but it can be easily changed.
Closest provides, in addition, ready to use, the metrics defined by the L1, L2,
Lp and LM norms.  New metric functions are very easy to write. In general,
closest methods work fine with distance function that are not *metrics* in the
mathematical sense, as long as they are positive definite: they need not to be
necessarily symmetric nor verify the triangle inequality. For example,
regarding a non-symmetric distance function, closest will always (and only)
call it with the query point as the first argument with the cell, brute-force
and culling methods. The kd-tree method will also call the distance function
while pruning the search tree, so results might differ from the other methods.
If in doubt, check the your solutions for a particular function comparing with
the brute-force method. 

For vector spaces with a norm, which method to use depends mainly on the number
of dimensions and the size of the problem. For low dimensions :math:`N_{D} < 6`
the fastest method is usually the cell method. If the number of points in the
set is large, :math:`N > 10^7` and :math:`N_{D} < 6`, or the number of
dimensions  :math:`6 < N_{D} \lesssim 20`, then the kd-tree method is the
generally the fastest. For a number of dimensions :math:`N_{D} > 20` the
brute-force method is generally faster than any other. The dimensional culling
method is the slowest, always, unless the point set is small and has a very
biased distribution which might throw off-track both the cell and kd-tree
methods.

The memory footprint, except for the dimensional culling method in which is
quite memory hungry, is roughly :math:`N_{\mathrm{data}}` doubles and
:math:`N_{\mathrm{data}}` integers, where :math:`N_{\mathrm{data}}` are the
:math:`N_{D}`-dimensional data against which the searches will be performed. 

Download
========
The latest Closest is version 2.0, and can be downloaded here: `closest-2.0.zip`_. 

The source repository is hosted at `github`_.

Licensing 
=========
Closest is released under the `GPL-3.0`_, but if you need it under other licensing 
terms, please feel free to `email me <mailto:trifling.github@gmail.com>`_.

Getting started
===============

Installing closest
------------------

Closest uses cmake, so it should be easy and fairly automatic to build. For example,
to build and install a dynamic version of closest in /usr/local do:

.. highlight:: bash

.. code-block:: bash

   $ mkdir bld 
   $ cd bld 
   $ cmake -DCMAKE_INSTALL_PREFIX=/usr/local ../ 
   $ make 
   $ make install 

Using Closest
-------------

Closest uses one joint C/C++ header file, closest.h, so it’s enough to put the line

.. code-block:: c

   #include <closest.h>

at the beginning of every C or C++ source file that uses Closest. For the Python interface
it is enough to

.. code-block:: python

   import closest

For the C API there is also just two libraries to link with, libclosest itself and the C
math library `m`. Compile and link a program using closest as follows:

.. code-block:: bash

   cc -o prog prog.c -lclosest -lm

The C++ API it also needs to be linked with the libclosest_cxx library:

.. code-block:: bash

   c++ -o prog prog.cpp -lclosest_cxx -lclosest -lm


Tutorials
=========

Simple example
--------------

Suppose we want to find the k=5 nearest neighbors to a point x in an array data 
of 10000 points in d=100 dimensions. Using the cell-based method:

First include the required header

.. code-block:: c

   #include<closest.h>   

Then make the appropriate calls

.. code-block:: c

   double *data;  /* points */
   int k = 5;     /* number of nearest neighbors to calculate */
   int d = 100;   /* number of dimensions */
   int n = 10000; /* number of data points */

   /* 
      assume data has been alloc'd filled with 10000 
      points in d=100 dimensions, each of the 100 
      coordinates for a point followed by the next 
      until all 10000 are specified.
      
      x is the point where the nearest neighbors 
      query is done.
   */
      
      
   /* initialize the search */
   cell_t *cell = cell_init( d, n, data, -1 );
   
   /* alloc space for the answer */
   int i_cell[k];
   int d_cell[k];

   /* 
      cell_knearest returns 
        in i_cell, the index for the points in the 
                   data array of the k-nearest neighbors to x 
        in d_cell, the distances for those same points to x 
   */
   cell_knearest( cell, x, k, i_cell, d_cell ); 

   /* 
      any other number of queries for different x could be done at 
      this point but we are done, so  free resources 
   */
   cell_free(cell);


Writing you own distance functions
----------------------------------

Writing custom distance functions is very easy. For example, a metric
based on the infinity norm

.. code-block:: c

   double inf( int nd, double *x, double *y, void *data ) {
      double t = DBL_MIN;
      for( int k=0; k<nd; k++ ) {
         double a = fabs(x[k]-y[k]);
         if( a > t )
            t = a;
      }
      return t;
   }

can be plugged in a breeze:

.. code-block:: c
   
   /* init */
   cell_t *cell = cell_init( d, n, data );
   
   /* set the custom metric */
   cell_set_metric( cell, inf, NULL );   

   /* search */
   cell_knearest( cell, x, k, i_cell, d_cell ); 

   /* clean up */
   cell_free(cell);

As stated in the introduction, Closest is not very picky about the
kind of distance functions one can plug and get correct results. However,
if in doubt, please check the results against the brute-force method
solutions. 

API Reference 
=============

C API
-----

There are four different sets of functions, one for each implemented algorithm,
cell, tree, cull-based and brute-force, with a similar interface.  To use them
`closest.h` needs to be included:

.. code-block:: c

   #include<closest.h>   

.. c:function:: cell_t *cell_init( int nd, int ni, double *xi, int ncpd )
                brut_t *brut_init( int nd, int ni, double *xi )
                cull_t *brut_init( int nd, int ni, double *xi, int prob )
                tree_t *brut_init( int nd, int ni, double *xi, int ppn )

   + int `nd`, input: the number of dimensions
   + int `ni`, input: the number of points
   + double `*xi`, input: pointer to the points, a \[`ni` x `nd`\] array in row-major order.
   + int `ncpd`, input: for the cell based method, the number of subdivision cells per dimension, 
     or -1 to let closest decide a reasonable value. Performance can potentially be improved by
     carefully choosing `ncpd`.
   + double `prob`, input: for the cull based method, a probability that defines the initial
     search radius of the algorithm. 0 < prob < 1. Pass -1.0 to take the default value. Performance 
     can potentially be improved by carefully choosing `prob`.
   + int `ppn`, input: what is generally know as bucket size, or the max number of points inside 
     each leaf-node in the tree decomposition. By default 32. Pass -1 to let closest decide. Performance 
     can potentially be improved by carefully choosing `ppn`.
   + returns: a cell_t, tree_t, cull_t or brut_t pointer

   Process the `ni` points `xi` in `nd`-dimensions such that any number of 
   k-nearest-neighbor searches can be performed subsequently against those points. 
   Returns a pointer to a data struct type that needs to be passed to \*_knearest set
   of functions. Once you are done searching, free up memory by calling 
   \*_free on the pointer.

.. c:function:: int cell_knearest( cell_t *t, int nq, double *xq, int no, int *index, double *distance )
                int brut_knearest( brut_t *t, int nq, double *xq, int no, int *index, double *distance )
                int cull_knearest( cull_t *t, int nq, double *xq, int no, int *index, double *distance )
                int tree_knearest( tree_t *t, int nq, double *xq, int no, int *index, double *distance )

   + cell_t, brut_t, cull_t or tree_t `*t`, input: the pointer returned by the \*_init set of functions
   + int `nq`, input: the number of query points
   + double `*xq`, input: a pointer to an [`nq` x `nd`] array with the coordinates of the query points 
   + int `no`, input: the number of nearest neighbors to find per query point
   + int `*index`: on input, a pointer to a [`nq` x `no`] array; on output it will contain the indices 
     of the `no` closest neighbors to the query points, ordered by distance for each.
   + double `*distance`: on input, a pointer to a [`nq` x `no`] array; on output
     it will contain the distances of the nearest neighbors found.
   + returns: the number of neighbors found, if successful equal to `nq` x `no`.
     
   Find the `no` nearest neighbors to the given query points `xq`, and return
   the indices to the original points array in `index` and the respective distances
   in `distance`.

.. c:function:: void cell_free( cell_t *t )
                void brut_free( brut_t *t )
                void cull_free( cull_t *t )
                void tree_free( tree_t *t )

   + cell_t, brut_t, cull_t or tree_t `*t`, input: the pointer returned by the \*_init functions
     
   Free the resources allocated by the \*_init functions. No further calls to 
   \*_knearest or \*_set_distfun can be performed.

.. c:function:: void cell_set_metric( cell_t *t, metric fun, void *data )
                void brut_set_metric( brut_t *t, metric fun, void *data )
                void cull_set_metric( cull_t *t, metric fun, void *data )
                void tree_set_metric( tree_t *t, metric fun, void *data )

   + cell_t, brut_t, cull_t or tree_t `*t`, input: the pointer returned by the \*_init functions
   + fun, input: a function with the following signature: double fun( int n, double \*x, double \*y, void \*data)
     that calculates the distance between two given points `x` and `y` of dimension `n`. 
   + data, intput: a pointer that will be passed unaltered each time fun is called and that
     can be used to pass extra data to the distance function.
     
.. c:function:: double L1metric( int n, double *x, double *y, void *data)

      Implements the L1, Manhattan or taxicab metric.

.. c:function:: double L2metric( int n, double *x, double *y, void *data)

      Implements the L2 or Euclidean metric.

.. c:function:: double Lpmetric( int n, double *x, double *y, void *data)

      Implements the Lp metric. It can be used also as a fractional semimetric with 0 < p < 1.
      The `p` parameter shall be passed by the void* data pointer.

.. c:function:: double LMmetric( int n, double *x, double *y, void *data)

      Implements the LM metric d(x,y)=x M y where M is a nd x nd matrix that
      has to be passed through the void* data pointer.

C++ API
-------

It is a very thin wrapper to the C API, which makes it ideal
to plug in conjunction with you preferred array classes as
it is framework-agnostic.

To use them `closest.h` needs to be included from a cpp file.

.. code-block:: c

   #include<closest.h>   

All is defined within the `closest` namespace. The metric functions are the
same as in the C interface, and it offers the following classes with the same
parameters as the C interface:


.. cpp:class:: closest::Cell 

.. cpp:function:: closest::Cell::Cell( int nd, int ni, double *xi, int npcd )

.. cpp:function:: closest::Cell::knearest( int nq, double *xq, int no, int *idx, double *dst )

.. cpp:function:: closest::Cell::set_metric( metric fun, void *data )


.. cpp:class:: closest::Brut

.. cpp:function:: closest::Brut::Brut( int nd, int ni, double *xi )

.. cpp:function:: closest::Brut::knearest( int nq, double *xq, int no, int *idx, double *dst )

.. cpp:function:: closest::Brut::set_metric( metric fun, void *data )


.. cpp:class:: closest::Cull

.. cpp:function:: closest::Cull::Cull( int nd, int ni, double *xi, double prob ) 

.. cpp:function:: closest::Cull::knearest( int nq, double *xq, int no, int *idx, double *dst )

.. cpp:function:: closest::Cull::set_metric( metric fun, void *data )


.. cpp:class:: closest::Tree 

.. cpp:function:: closest::Tree::Tree( int nd, int ni, double *xi, int npcd )

.. cpp:function:: closest::Tree::knearest( int nq, double *xq, int no, int *idx, double *dst )

.. cpp:function:: closest::Tree::set_metric( metric fun, void *data )



Python API
----------

To use it the `closest` module needs to be imported:

.. code-block:: python

   import closest   

.. py:class:: Cell( data, ncpd = -1 ) 
   
   + `data`, input: a numpy array with shape (ni,nd), where ni is the number
     of data points and nd is the number of dimensions.
   + int `ncpd`, input: the number of subdivision cells per dimension, or -1 to 
     let closest decide a reasonable value. Performance can potentially be improved by
     carefully choosing ncpd.

   Process the points in data such that any number of k-nearest-neighbor searches can 
   be performed subsequently against those points. 

.. py:method:: Cell.knearest( x, n, metric = None, metric_data = None )

   + `x`, input: a numpy array with shape (nq,nd), where nq is the number
     of different query points and nd in the number of dimensions.
   + metric: None (L2 will be used), or a string with 'L1', 'L2' or 'Lp' to use those 
     metrics, or a python function that returns a double and has the signature 
     py_metric( `x`, `y`, `data` ), where `x` and `y` are numpy arrays of shape (`nd`,) being `nd` the number of dimensions, and `data`
     a python object with extra data or None.  
   + metric_data: any python object that will be passed unaltered to the metric function 

   Find the `n` nearest neighbors to the given query points `x`, and return
   a tuple (idx,dst) with the indices to the original points array in `idx` 
   and the respective distances in `dst`.

The brute-force, cull and tree interfaces are very similar, with the following signatures:

.. py:class:: Brut( data ) 

.. py:method:: Brut.knearest( x, n, metric = None, metric_data = None )

.. py:class:: Cull( data, prob = -1.0 ) 

.. py:method:: Cull.knearest( x, n, metric = None, metric_data = None )

.. py:class:: Tree( data, ppn = -1 ) 

.. py:method:: Tree.knearest( x, n, metric = None, metric_data = None )


