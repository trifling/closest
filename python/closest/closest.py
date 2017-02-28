#  
#  This file is part of closest, a library for cull-based and 
#  cell-based spatial search of k-nearest neighbors in n-dimensions.
#  Copyright (C) 2011-2016 Daniel Pena <trifling.github@gmail.com>
#
#  Closest is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Closest is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with closest. If not, see <http://www.gnu.org/licenses/>.
# 

import numpy as npy
import os

from ctypes import CDLL, POINTER, c_int, c_double, byref, c_void_p, Structure, util
from ctypes.util import find_library

found = False
lpath = find_library( 'closest' )
if not lpath:
    try:
        lpath = os.path.join( os.environ['CLOSEST_PATH'], 'libclosest.so' )
    except:
        pass
    try:
        for root in os.environ['LD_LIBRARY_PATH'].split(':'):
            lpath = os.path.join( root, 'libclosest.so' )
            if os.path.isfile(lpath):
                break
    except:
        pass

if not lpath:
    raise RuntimeError( 'could not locate the closest library, please set CLOSEST_PATH to the valid location')

_closest = CDLL( lpath )

_cell_init = _closest.cell_init
_cell_init.restype = c_void_p
_cell_knearest = _closest.cell_knearest
_cell_knearest.restype = c_int
_cell_free = _closest.cell_free
_cell_free.restype = None

_cull_init = _closest.cull_init
_cull_init.restype = c_void_p
_cull_knearest = _closest.cull_knearest
_cull_knearest.restype = c_int
_cull_free = _closest.cull_free
_cull_free.restype = None

_bruteforce_knearest = _closest.bruteforce_knearest
_bruteforce_knearest.restype = c_int 

class Cell:
    def __init__( self, data ):
        if data.dtype != 'float64':
            raise TypeError( 'closest Cell class takes only float64 numpy arrays' )

        self.data = data
        self.ni = data.shape[0]
        self.nd = data.shape[1]

        nr = c_int(-1)
        nd = c_int(self.nd)
        ni = c_int(self.ni)
        pr = data.ctypes.data_as( POINTER(c_double) )
        self.handle = _cell_init( nd, ni, nr, pr )

    def __del__( self ):
        i = _cell_free( self.handle )

    def knearest( self, x, n ):
        no = c_int(n)
        pr = x.ctypes.data_as( POINTER(c_double) )
        
        idx = npy.zeros( n, dtype=npy.int32 )
        dst = npy.zeros( n, dtype=npy.float64 )
        
        pidx = idx.ctypes.data_as(POINTER(c_int))
        pdst = dst.ctypes.data_as(POINTER(c_double))

        ret = _cell_knearest( self.handle, pr, no, pidx, pdst )

        return ret, idx, dst

class Cull:
    def __init__( self, data ):
        if data.dtype != 'float64':
            raise TypeError( 'closest Cull class takes only float64 numpy arrays' )

        self.data = data
        self.ni = data.shape[0]
        self.nd = data.shape[1]

        nd = c_int(self.nd)
        ni = c_int(self.ni)
        pr = data.ctypes.data_as( POINTER(c_double) )
        self.handle = _cull_init( nd, ni, pr )

    def __del__( self ):
        i = _cull_free( self.handle )

    def knearest( self, x, n ):
        no = c_int(n)
        pr = x.ctypes.data_as( POINTER(c_double) )
        
        idx = npy.zeros( n, dtype=npy.int32 )
        dst = npy.zeros( n, dtype=npy.float64 )
        
        pidx = idx.ctypes.data_as(POINTER(c_int))
        pdst = dst.ctypes.data_as(POINTER(c_double))

        ret = _cull_knearest( self.handle, pr, no, pidx, pdst )

        return ret, idx, dst

class Bruteforce:
    def __init__( self, data ):
        if data.dtype != 'float64':
            raise TypeError( 'closest Bruteforce class takes only float64 numpy arrays' )

        self.data = data
        self.ni = data.shape[0]
        self.nd = data.shape[1]

    def knearest( self, x, n ):

        nd = c_int(self.nd)
        ni = c_int(self.ni)
        xi = self.data.ctypes.data_as( POINTER(c_double) )

        no = c_int(n)
        pr = x.ctypes.data_as( POINTER(c_double) )
        
        idx = npy.zeros( n, dtype=npy.int32 )
        dst = npy.zeros( n, dtype=npy.float64 )
        
        pidx = idx.ctypes.data_as(POINTER(c_int))
        pdst = dst.ctypes.data_as(POINTER(c_double))

        ret = _bruteforce_knearest( nd, ni, xi, pr, no, pidx, pdst )

        return ret, idx, dst

if __name__ == '__main__':
    
    nd = 5
    ni = 20
    data = npy.random.rand( ni, nd )
    x = npy.random.rand( nd )
    print( 'DATA=')
    print( data )
    cell = Cell(data)
    print()
    print( '  X=',x )
    print()
    ret, idx, dst = cell.knearest( x, 3 )
    print( 'RET=', ret )
    print( 'IDX=', idx )
    print( 'DST=', dst )
    print( 'LST=' )
    for i in idx:
        print( data[i,:] )

    print()
    bfor = Bruteforce(data)
    ret, idx, dst = bfor.knearest( x, 3 )
    for i in idx:
        print( data[i,:] )

