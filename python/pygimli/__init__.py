# -*- coding: utf-8 -*-

"""
These are the python bindings for libgimli

import pygimli

or

import pygimli as g

or

from pygimli import *
"""

import sys,os
if sys.platform == 'win32':
    os.environ['PATH'] =  __path__[0] +';' + os.environ['PATH']

try:
    from _pygimli_ import *
except ImportError as e:
    print e
    sys.stderr.write("ERROR: cannot import the library '_pygimli_'.\n")

import locale
if locale.localeconv()['decimal_point'] == ',':
    print "Found locale decimal_point ',', change it to: decimal point '.':",
    locale.localeconv()['decimal_point']
    locale.setlocale( locale.LC_NUMERIC, 'C' )

############################
# print function for gimli stuff
############################
def RVector_str( self ):
    s = str( type( self ) ) + " "  + str( self.size() );

    if len( self ) == 0:
        return s
    else:
        s += " ["

    if len( self ) < 51:
        for i in range( 0, len( self )-1 ):
            s = s + str( self[ i ] ) + ", "

        s = s + str( self[ len( self )-1 ] ) + "]"
        return s
    return str( type( self ) ) + " " + str( self.size() ) + " [" + str( self[ 0 ] ) + ",...," + str( self[ self.size() -1 ] ) + "]"


def RVector3_str( self ):
    return "RVector3: (" + str( self.x() ) + ", " + str( self.y() ) + ", " + str( self.z() ) + ")"

def RMatrix_str( self ):
    return "RMatrix: " + str( self.rows() ) + " x " + str( self.cols() )

def Line_str( self ):
    return "Line: " + str( self.p0() ) + "  " + str( self.p1() )

def Mesh_str( self ):
    return "Mesh: Nodes: " + str( self.nodeCount() ) + " Cells: " + str( self.cellCount() ) + " Boundaries: " + str( self.boundaryCount() )

def Data_str( self ):
    return "Data: Sensors: " + str( self.sensorCount() ) + " data: " + str( self.size() )

_pygimli_.RVector3.__str__ = RVector3_str
_pygimli_.RVector.__str__ = RVector_str
_pygimli_.BVector.__str__ = RVector_str
_pygimli_.RMatrix.__str__ = RMatrix_str
_pygimli_.Line.__str__ = Line_str
_pygimli_.Mesh.__str__ = Mesh_str
_pygimli_.DataContainer.__str__ = Data_str
_pygimli_.stdVectorUL.size = _pygimli_.stdVectorUL.__len__
_pygimli_.stdVectorUL.__str__ = RVector_str

############################
# compatibility stuff
############################
def nonzero_test( self ):
    raise "Warning! there is no 'and' for BVector and RVector use '&' instead"

_pygimli_.RVector.__nonzero__ = nonzero_test
_pygimli_.BVector.__nonzero__ = nonzero_test

############################
# allow:
############################
def __ADD( self, val ):
    ret = type( self )()
    for i, r in enumerate( self ): ret.append( r + val )
    return ret

_pygimli_.stdVectorUL.__add__ = __ADD

############################
# Indexing operator for RVector, RVector3, RMatrix
############################
def __getVal( self, idx ):
    if type( idx ) is slice:
        if idx.step is None:
            return self( long( idx.start ), long( idx.stop ) )
        else:
            "not yet implemented"

    if idx == -1: idx = len( self )-1
    
    return self.getVal( long( idx ) )
# def __getVal( ... )
    
def __setVal( self, idx, val ):
    self.setVal( val, idx )

_pygimli_.RVector.__setitem__ = __setVal
_pygimli_.RVector.__getitem__ = __getVal

_pygimli_.RVector3.__setitem__ = __setVal

def __getValMatrix( self, idx ):
    if idx == -1: idx = len( self )-1
    return self.rowR( idx )

_pygimli_.RMatrix.__getitem__ = __getValMatrix
_pygimli_.RMatrix.__setitem__ = __setVal


############################
# len( RVector ), RMatrix
############################
def RVector_len( self ): return self.size()

_pygimli_.RVector.__len__ = RVector_len
_pygimli_.BVector.__len__ = RVector_len

def RMatrix_len( self ): return self.rows()
_pygimli_.RMatrix.__len__ = RMatrix_len


############################
# Iterator support for RVector allow to apply python build-ins
############################
class VectorIter:
    def __init__( self, vec ):
        self.vec = vec
        self.length = len( vec )
        self.pos = -1;

    def __iter__( self ): return self

    def next( self ):
        self.pos += 1
        if self.pos == self.length:
            raise StopIteration()
        else:
            return self.vec[ self.pos ]

def __VectorIterCall__( self ): return VectorIter( self )
_pygimli_.RVector.__iter__  = __VectorIterCall__
_pygimli_.BVector.__iter__  = __VectorIterCall__
_pygimli_.RVector3.__iter__ = __VectorIterCall__
_pygimli_.RMatrix.__iter__  = __VectorIterCall__

############################
# non automatic exposed functoions
############################
def abs( v ):
    return fabs( v )


############################
# ???
############################
def RVector_ArrayInit( self ):
    #self.ndim = [ 1, self.size() ]
    self.ndim = 1
_pygimli_.RVector.isArray = RVector_ArrayInit


#
# Construct RVector from numpy array , opposite to asarray( RVector ) 
#
def asvector( arr ):
    r = _pygimli_.RVector( len( arr ) )
    
    for i, v in enumerate( arr ):
        r[ i ] = v

    return r
    

###########################
# unsorted stuff
###########################

#def __interpolate__cmd( self ):
    #print "my interpolate"
    #return _pygimli_.interpolate_GILsave__
#_pygimli_.interpolate = __interpolate__cmd

_pygimli_.interpolate = _pygimli_.interpolate_GILsave__

# define template void Quaternion< double >::rotMatrix( Matrix < double > & mat ) const;
def __getRotMatrix__( self, mat) : getRotMatrix__( self, mat )
_pygimli_.RQuaternion.rotMatrix = __getRotMatrix__


#some rvector helpers
def randN( self, n ):
    '''Create RVector of length n with normally distributed random numbers'''
    r = _pygimli_.RVector( n )
    _pygimli_.randn( r )
    return r

def checkAndFixLocaleDecimal_point( verbose = False ):
    import locale
    if locale.localeconv()['decimal_point'] == ',':
        if verbose:
            print "decimal point: ", locale.localeconv()['decimal_point']
            print "setting ."
        locale.setlocale( locale.LC_NUMERIC, 'C' )
