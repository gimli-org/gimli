# -*- coding: utf-8 -*-

"""
These are the python bindings for libgimli

import pygimli

or

import pygimli as pg

or

from pygimli import *
"""
from __future__ import print_function

import os
import subprocess
import sys

if sys.platform == 'win32':
    os.environ['PATH'] = __path__[0] + ';' + os.environ['PATH']

try:
    from . _pygimli_ import *
except ImportError as e:
    print(e)
    import traceback

    traceback.print_exc(file=sys.stdout)
    sys.stderr.write("ERROR: cannot import the library '_pygimli_'.\n")

import locale
if locale.localeconv()['decimal_point'] == ',':
    print(
        "Found locale decimal_point ',', change it to: decimal point '.':",
        end=' ')
    locale.localeconv()['decimal_point']
    locale.setlocale(locale.LC_NUMERIC, 'C')

###  Global shortcutes #####
_pygimli_.load = None
from pygimli.io import load
from pygimli.viewer import show
############################


############################
# print function for gimli stuff
############################
def RVector_str(self, valsOnly=False):
    s = str()

    if not valsOnly:
        s = str(type(self)) + " " + str(self.size())

    if len(self) == 0:
        return s
    else:
        s += " ["

    if len(self) < 101:
        for i in range(0, len(self) - 1):
            s = s + str(self[i]) + ", "

        s = s + str(self[len(self) - 1]) + "]"
        return s
    return (
        str(type(self)) + " " + str(self.size()) +
        " [" + str(self[0]) + ",...," + str(self[self.size() - 1]) + "]"
    )


def RVector3_str(self):
    return (
        "RVector3: (" + str(self.x()) + ", " +
        str(self.y()) + ", " + str(self.z()) + ")"
    )


def RMatrix_str(self):
    s = "RMatrix: " + str(self.rows()) + " x " + str(self.cols())

    if (self.rows() < 6):
        s += '\n'
        for v in range(self.rows()):
            s += self[v].__str__(True) + '\n'
    return s


def Line_str(self):
    return "Line: " + str(self.p0()) + "  " + str(self.p1())


def Mesh_str(self):
    return (
        "Mesh: Nodes: " + str(self.nodeCount()) + " Cells: " +
        str(self.cellCount()) + " Boundaries: " +
        str(self.boundaryCount())
    )


def Data_str(self):
    return (
        "Data: Sensors: " +
        str(self.sensorCount()) + " data: " + str(self.size())
    )

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

def nonzero_test(self):
    raise BaseException("Warning! there is no 'and' and 'or' for BVector and RVector. " + \
        "Use binary operators '&' or '|' instead. " + \
        "If you looking for the nonzero test, use len(v) > 0")

_pygimli_.RVector.__nonzero__ = nonzero_test
_pygimli_.BVector.__nonzero__ = nonzero_test
_pygimli_.RVector.__bool__ = nonzero_test
_pygimli_.BVector.__bool__ = nonzero_test


############################
# allow:
############################

def __ADD(self, val):
    ret = type(self)()
    for i, r in enumerate(self):
        ret.append(r + val)
    return ret

_pygimli_.stdVectorUL.__add__ = __ADD


############################
# Indexing operator for RVector, RVector3, RMatrix
############################
def __getVal(self, idx):
    """
        Hell slow
    """
    
    if isinstance(idx, list):
        idxL = _pygimli_.stdVectorUL()
        for ix in idx:
            idxL.append(ix)
        return self(idxL)

    if isinstance(idx, stdVectorUL):
        return self(idx)

    if isinstance(idx, slice):
        if idx.step is None:
            return self(int(idx.start), int(idx.stop))
        else:
            ids = range(idx.start, idx.stop, idx.step)
            if len(ids):
                return self(ids)
            else:
                raise Exception("slice invalid")

    if idx == -1:
        idx = len(self) - 1

    return self.getVal(int(idx))
# def __getVal(...)


def __setVal(self, idx, val):

    if isinstance(idx, slice):
        if idx.step is None:
            self.setVal(float(val), int(idx.start), int(idx.stop))
            return
        else:
            "not yet implemented"

    self.setVal(val, idx)

def __getValMatrix(self, idx):

    if isinstance(idx, slice):
        if idx.step is None:

            start = idx.start
            if idx.start is None:
                start = 0

            stop = idx.stop
            if idx.stop is None:
                stop = len(self)

            return [self.rowR(i) for i in range(start, stop)]

            # return self(long(idx.start), long(idx.stop))
        else:
            "not yet implemented"

    if idx == -1:
        idx = len(self) - 1

    return self.rowR(idx)


_pygimli_.RVector.__setitem__ = __setVal
_pygimli_.RVector.__getitem__ = __getVal # very slow -- inline is better

_pygimli_.RVector3.__setitem__ = __setVal

_pygimli_.RMatrix.__getitem__ = __getValMatrix # very slow -- inline is better
_pygimli_.RMatrix.__setitem__ = __setVal


############################
# len(RVector), RMatrix
############################
def RVector_len(self):
    return self.size()

_pygimli_.RVector.__len__ = RVector_len
_pygimli_.BVector.__len__ = RVector_len


def RMatrix_len(self):
    return self.rows()
_pygimli_.RMatrix.__len__ = RMatrix_len


############################
# Iterator support for RVector allow to apply python build-ins
############################
class VectorIter:

    def __init__(self, vec):
        self.vec = vec
        self.length = len(vec)
        self.pos = -1

    def __iter__(self):
        return self

    # this is for python < 3
    def next(self):
        return self.__next__()

    # this is the same but for python > 3
    def __next__(self):
        #print('next')
        self.pos += 1
        if self.pos == self.length:
            raise StopIteration()
        else:
            return self.vec[self.pos]

def __VectorIterCall__(self):
    #return VectorIter(self)
    return _pygimli_.RVectorIter(self.beginPyIter())

def __MatIterCall__(self):
    return VectorIter(self)

class Vector3Iter (VectorIter):

    def __init__(self, vec):
        self.vec = vec
        self.length = 3
        self.pos = -1

def __Vector3IterCall__(self):
    return Vector3Iter(self)

_pygimli_.RVector.__iter__ = __VectorIterCall__
_pygimli_.BVector.__iter__ = __VectorIterCall__
_pygimli_.RVector3.__iter__ = __Vector3IterCall__
_pygimli_.RMatrix.__iter__ = __MatIterCall__

# DEPRECATED for backward compatibility should be removed
def asvector(array):
    return array
 
########## c to python converter ######
# default converter from RVector3 to numpy array 
def __RVector3ArrayCall__(self, idx=None):
    if idx:
        print(self)
        print(idx)
        raise Exception("we need to fix this")
    import numpy as np
    return np.array([self.getVal(0), self.getVal(1), self.getVal(2)])

# default converter from RVector to numpy array 
def __RVectorArrayCall__(self, idx=None):
    if idx:
        print(self)
        print(idx)
        raise Exception("we need to fix this")
    return self.array()
    
_pygimli_.RVector.__array__ = __RVectorArrayCall__
_pygimli_.RVector3.__array__ = __RVector3ArrayCall__
#_pygimli_.RVector3.__array__ = _pygimli_.RVector3.array
#del _pygimli_.RVector.__array__

############################
# non automatic exposed functions
############################
def abs(v):
    return fabs(v)


########################################################
# compare operators for stdVector
########################################################

def __CMP_stdVectorI__(self, val):
    raise Exception("__CMP_stdVectorI__ do not use")
    return 0
    print("__CMP_stdVectorI__")
    ret = _pygimli_.BVector(len(self))
    for i, v in enumerate(ret):
        print(self[i] < val, int(self[i] < val))
        #v = int(self[i] < val)
    print(ret)


def __EQ_stdVectorI__(self, val):
    raise Exception("__EQ_stdVectorI__ do not use")
    ret = _pygimli_.BVector(len(self))
    for i, v in enumerate(ret):
        print(self[i] == val, int(self[i] == val))
        #v = self[i] == val
    print(ret)

_pygimli_.stdVectorI.__cmp__ = __CMP_stdVectorI__
_pygimli_.stdVectorI.__eq__ = __EQ_stdVectorI__


############################
# usefull stuff
############################


############################
# ???
############################
#def RVector_ArrayInit(self):
    ##self.ndim = [1, self.size()]
    #self.ndim = 1
#_pygimli_.RVector.isArray = RVector_ArrayInit


#
# Construct RVector from numpy array , opposite to asarray(RVector)
#
#def asvector(arr):
    #r = _pygimli_.RVector(len(arr))

    #for i, v in enumerate(arr):
        #r[i] = v

    #return r

# PEP conform version str with SVN revision number
def __svnversion__(path=__path__[0]):
    p = subprocess.Popen("svnversion -n %s" % path, shell=True, \
       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdout, stderr) = p.communicate()
    return str(stdout)

__version__ = _pygimli_.versionStr() + "_rev" + __svnversion__()

###########################
# unsorted stuff
###########################

# def __interpolate__cmd(self):
    # print "my interpolate"
    # return _pygimli_.interpolate_GILsave__
#_pygimli_.interpolate = __interpolate__cmd

_pygimli_.interpolate = _pygimli_.interpolate_GILsave__

# define template void Quaternion< double >::rotMatrix(Matrix < double > &
# mat) const;


def __getRotMatrix__(self, mat):
    getRotMatrix__(self, mat)
_pygimli_.RQuaternion.rotMatrix = __getRotMatrix__


# some rvector helpers
def randN(self, n):
    '''Create RVector of length n with normally distributed random numbers'''
    r = _pygimli_.RVector(n)
    _pygimli_.randn(r)
    return r


def checkAndFixLocaleDecimal_point(verbose=False):
    import locale
    if locale.localeconv()['decimal_point'] == ',':
        if verbose:
            print("decimal point: ", locale.localeconv()['decimal_point'])
            print("setting .")
        locale.setlocale(locale.LC_NUMERIC, 'C')
