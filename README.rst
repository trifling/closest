
.. image:: https://travis-ci.org/trifling/closest.svg?branch=master
   :target: https://travis-ci.org/trifling/closest

.. _`is well documented`: http://trifling-matters.com/closest.html
.. _`GPL-3.0`: https://opensource.org/licenses/GPL-3.0

Closest README
==============

A C library for k-nearest neighbor search in N-dimensions with arbitrary
metric functions. 

Closest:

- is simple to use
- is flexible (custom distance functions)
- provides API interfaces for C, C++ and Python
- written in ISO C99
- has no external dependencies
- is small and fast

- `is well documented`_

Closest is licensed under the `GPL-3.0`_ (see the LICENSE file for details).

Code Example
============

Find the `k` = 10 nearest neighbors to a point `x` in an array `data` of 1000 
points in `d` = 20 dimensions:

.. code-block:: c

   int d = 20;
   int k = 10;
   int n = 1000;

   /* 
      double *data[n*d] has been filled with 1000 points in d=20 dimensions, each
      of the 20 coordinates for a point followed by the next until all 1000 
      are specified.
     
      x is the point where the nearest neighbors query is done
   */

   /* initialize the search */
   cell_t *cell = cell_init( d, n, data, -1 );
   
   int i_cell[k];
   int d_cell[k];
   /* cell_knearest returns 
        in i_cell, the index for the points in the `data` array of the k-nearest neighbors to x 
        in d_cell, the distances for those same points to x */
   cell_knearest( cell, 1, x, k, i_cell, d_cell ); 
   cell_free(cell);

   
Compilation and installation
============================

Closest uses cmake as build system. To build and install a dynamic version of closest in /usr/local:

.. code-block:: bash
   $ mkdir bld 
   $ cd bld 
   $ cmake -DCMAKE_INSTALL_PREFIX=/usr/local ../ 
   $ make 
   $ make install


Documentation
-------------

Documentation is available at http://trifling-matters.com/closest.html.

The documentation source is in the ``doc/`` subdirectory. To generate
the HTML documentation, invoke::

   $ make html

Then, point your browser to ``doc/bld/index.html``. 

