#  This file is part of closest, a library for k-nearest neighbors (kNN) search
#  in n-dimensions.  
#  Copyright (C) 2011-2016 Daniel Pena <trifling.github@gmail.com>
#
#  Closest is free software: you can redistribute it and/or modify it under the
#  terms of the GNU General Public License as published by the Free Software
#  Foundation, either version 3 of the License, or (at your option) any later
#  version.
#
#  Closest is distributed in the hope that it will be useful, but WITHOUT ANY
#  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
#  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
#  details.
#
#  You should have received a copy of the GNU General Public License along with
#  closest. If not, see <http://www.gnu.org/licenses/>.

import numpy as npy
import os

from ctypes import CDLL, POINTER, c_int, c_double, byref, c_void_p, Structure, util, CFUNCTYPE, py_object
from ctypes.util import find_library

_found = False
_lpath = find_library( 'closest' )
if not _lpath:
    try:
        _lpath = os.path.join( os.environ['CLOSEST_PATH'], 'libclosest.so' )
    except:
        pass
    try:
        for root in os.environ['LD_LIBRARY_PATH'].split(':'):
            _lpath = os.path.join( root, 'libclosest.so' )
            if os.path.isfile(_lpath):
                break
    except:
        pass

if not _lpath:
    raise RuntimeError( 'could not locate the closest library, please set CLOSEST_PATH to the valid location')

_closest = CDLL( _lpath )


_cull_init = _closest.cull_init
_cull_init.restype = c_void_p
_cull_knearest = _closest.cull_knearest
_cull_knearest.restype = c_int
_cull_free = _closest.cull_free
_cull_free.restype = None
_cull_set_metric = _closest.cull_set_metric
_cull_set_metric.restype = None

_cell_init = _closest.cell_init
_cell_init.restype = c_void_p
_cell_knearest = _closest.cell_knearest
_cell_knearest.restype = c_int
_cell_free = _closest.cell_free
_cell_free.restype = None
_cell_set_metric = _closest.cell_set_metric
_cell_set_metric.restype = None

_brut_init = _closest.brut_init
_brut_init.restype = c_void_p
_brut_knearest = _closest.brut_knearest
_brut_knearest.restype = c_int
_brut_free = _closest.brut_free
_brut_free.restype = None
_brut_set_metric = _closest.brut_set_metric
_brut_set_metric.restype = None

_tree_init = _closest.tree_init
_tree_init.restype = c_void_p
_tree_knearest = _closest.tree_knearest
_tree_knearest.restype = c_int
_tree_free = _closest.tree_free
_tree_free.restype = None
_tree_set_metric = _closest.tree_set_metric
_tree_set_metric.restype = None

_L1metric = _closest.L1metric
_L2metric = _closest.L2metric
_Lpmetric = _closest.Lpmetric
_LMmetric = _closest.LMmetric

def metric_wrap( py_metric, py_data ): 
    CDISTANCETYPE  = CFUNCTYPE( c_double, c_int, POINTER(c_double), POINTER(c_double), POINTER(py_object) ) 
    def c_metric_wrapper( n, px, py, pdata ):
        x = npy.ctypeslib.as_array(px,shape=(n,)) 
        y = npy.ctypeslib.as_array(py,shape=(n,)) 
        data = pdata.contents.value
        return py_metric( x, y, data )

    if py_metric is None:
        c_fun = None
        c_data = None

    elif type(py_metric)==str:
        if py_metric == 'L1':
            c_fun  = _L1metric
            c_data = None
        elif py_metric == 'L2':
            c_fun  = _L2metric
            c_data = None
        elif py_metric == 'Lp':
            c_fun  = _Lpmetric
            c_data = byref( c_double(py_data) )
        elif py_metric == 'LM':
            c_fun  = _LMmetric
            c_data = py_data.ctypes.data_as(POINTER(c_double))
        else:
            raise RuntimeError("Distance function " + py_metric + " not implemented")
    else:
        c_fun = CDISTANCETYPE(c_metric_wrapper)
        c_data = byref( py_object(py_data) )

    return c_fun, c_data

class Cull:
    def __init__( self, data, prob = -1 ):
        if data.dtype != 'float64':
            raise TypeError( 'closest Cull class takes only float64 numpy arrays' )

        self.data = data
        self.ni = data.shape[0]
        self.nd = data.shape[1]
        self.prob = prob

        nd = c_int(self.nd)
        ni = c_int(self.ni)
        pr = data.ctypes.data_as( POINTER(c_double) )
        ii = c_int(self.prob)
        self.handle = _cull_init( nd, ni, pr, ii )

    def __del__( self ):
        i = _cull_free( self.handle )

    def knearest( self, x, n, metric = None, distdata = None ):
        
        Lfun, Ldata = metric_wrap( metric, distdata )
        _cull_set_metric( self.handle, Lfun, Ldata )

        xx = x.flatten() # force numpy to put it in memory in (nd,nq) order
        nq = c_int( int( xx.size/self.nd ) )
        xq = xx.ctypes.data_as( POINTER(c_double) )
        no = c_int(n)
        idx = npy.zeros( (nq.value,no.value), dtype=npy.int32   )
        dst = npy.zeros( (nq.value,no.value), dtype=npy.float64 )
        
        pidx = idx.ctypes.data_as(POINTER(c_int))
        pdst = dst.ctypes.data_as(POINTER(c_double))

        ret = _cull_knearest( self.handle, nq, xq, no, pidx, pdst )
        
        return idx, dst 

class Cell:
    def __init__( self, data, ncpd = -1 ):
        if data.dtype != 'float64':
            raise TypeError( 'closest Cell class takes only float64 numpy arrays' )

        self.data = data
        self.ni = data.shape[0]
        self.nd = data.shape[1]
        self.ncpd = ncpd

        nd = c_int(self.nd)
        ni = c_int(self.ni)
        pr = data.ctypes.data_as( POINTER(c_double) )
        ii = c_int(self.ncpd)
        self.handle = _cell_init( nd, ni, pr, ii )

    def __del__( self ):
        i = _cell_free( self.handle )

    def knearest( self, x, n, metric = None, distdata = None ):
        
        Lfun, Ldata = metric_wrap( metric, distdata )
        _cell_set_metric( self.handle, Lfun, Ldata )

        xx = x.flatten() # force numpy to put it in memory in (nd,nq) order
        nq = c_int( int( xx.size/self.nd ) )
        xq = xx.ctypes.data_as( POINTER(c_double) )
        no = c_int(n)
        idx = npy.zeros( (nq.value,no.value), dtype=npy.int32   )
        dst = npy.zeros( (nq.value,no.value), dtype=npy.float64 )
        
        pidx = idx.ctypes.data_as(POINTER(c_int))
        pdst = dst.ctypes.data_as(POINTER(c_double))

        ret = _cell_knearest( self.handle, nq, xq, no, pidx, pdst )
        
        return idx, dst 

class Brut:
    def __init__( self, data ):
        if data.dtype != 'float64':
            raise TypeError( 'closest Bruteforce class takes only float64 numpy arrays' )

        self.data = data
        self.ni = data.shape[0]
        self.nd = data.shape[1]

        nd = c_int(self.nd)
        ni = c_int(self.ni)
        pr = data.ctypes.data_as( POINTER(c_double) )
        self.handle = _brut_init( nd, ni, pr )

    def __del__( self ):
        i = _brut_free( self.handle )

    def knearest( self, x, n, metric = None, distdata = None ):

        Lfun, Ldata = metric_wrap( metric, distdata )
        _brut_set_metric( self.handle, Lfun, Ldata )

        xx = x.flatten() # force numpy to put it in memory in (nd,nq) order
        nq = c_int( int( xx.size/self.nd ) )
        xq = xx.ctypes.data_as( POINTER(c_double) )
        no = c_int(n)
        idx = npy.zeros( (nq.value,no.value), dtype=npy.int32   )
        dst = npy.zeros( (nq.value,no.value), dtype=npy.float64 )
 
        pidx = idx.ctypes.data_as(POINTER(c_int))
        pdst = dst.ctypes.data_as(POINTER(c_double))

        ret = _brut_knearest( self.handle, nq, xq, no, pidx, pdst )

        return idx, dst 

class Tree:
    def __init__( self, data, ppn = -1 ):
        if data.dtype != 'float64':
            raise TypeError( 'closest Tree class takes only float64 numpy arrays' )

        self.data = data
        self.ni = data.shape[0]
        self.nd = data.shape[1]
        self.ppn = ppn

        nd = c_int(self.nd)
        ni = c_int(self.ni)
        pr = data.ctypes.data_as( POINTER(c_double) )
        ii = c_int(self.ppn)
        self.handle = _tree_init( nd, ni, pr, ii )

    def __del__( self ):
        i = _tree_free( self.handle )

    def knearest( self, x, n, metric = None, distdata = None ):
        
        Lfun, Ldata = metric_wrap( metric, distdata )
        _tree_set_metric( self.handle, Lfun, Ldata )

        xx = x.flatten() # force numpy to put it in memory in (nd,nq) order
        nq = c_int( int( xx.size/self.nd ) )
        xq = xx.ctypes.data_as( POINTER(c_double) )
        no = c_int(n)
        idx = npy.zeros( (nq.value,no.value), dtype=npy.int32   )
        dst = npy.zeros( (nq.value,no.value), dtype=npy.float64 )
        
        pidx = idx.ctypes.data_as(POINTER(c_int))
        pdst = dst.ctypes.data_as(POINTER(c_double))

        ret = _tree_knearest( self.handle, nq, xq, no, pidx, pdst )
        
        return idx, dst 

