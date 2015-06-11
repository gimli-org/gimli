# -*- coding: utf-8 -*-

"""
pyGIMLi - Python package for GIMLi including bindings to C++ library.

Usage:

.. code-block:: python

    import pygimli as pg

"""

from __future__ import print_function

import os
try:
    import subprocess  # check for 3.4
except ImportError:
    pass

import sys

if sys.platform == 'win32':
    os.environ['PATH'] = __path__[0] + ';' + os.environ['PATH']

_pygimli_ = None

try:
    from . import _pygimli_
    from . _pygimli_ import *
except ImportError as e:
    print(e)
    import traceback

    traceback.print_exc(file=sys.stdout)
    sys.stderr.write("ERROR: cannot import the library '_pygimli_'.\n")

import locale


def checkAndFixLocaleDecimal_point(verbose=False):
    if locale.localeconv()['decimal_point'] == ',':
        if verbose:
            print("Found locale decimal_point ',' "
                  "and change it to: decimal point '.'")
    try:
        locale.localeconv()['decimal_point']
        locale.setlocale(locale.LC_NUMERIC, 'C')
    except Exception as e:
        print(e)
        print('cannot set locale to decimal point')

checkAndFixLocaleDecimal_point(verbose=True)
# print(locale.localeconv()['decimal_point'])
# if locale.localeconv()['decimal_point'] == ',':
#    print("Found locale decimal_point ',' and change it to: decimal point '.'")
# try:
#    locale.localeconv()['decimal_point']
#    locale.setlocale(locale.LC_NUMERIC, 'C')
# except:
#    print('cannot set locale to decimal point')

###########################
###  Global shortcuts #####
###########################

_pygimli_.load = None
from pygimli.io import load
from pygimli.viewer import show, showLater, plt, wait
from pygimli.solver import solve

def showNow():
    showLater(0)

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
    return ("RVector3: (" + str(self.x()) + ", " +
            str(self.y()) + ", " + str(self.z()) + ")" )

def R3Vector_str(self):
    return ("R3Vector: n=" + str(self.size()))

def RMatrix_str(self):
    s = "RMatrix: " + str(self.rows()) + " x " + str(self.cols())

    if (self.rows() < 6):
        s += '\n'
        for v in range(self.rows()):
            s += self[v].__str__(True) + '\n'
    return s

def CMatrix_str(self):
    s = "CMatrix: " + str(self.rows()) + " x " + str(self.cols())

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

_pygimli_.RVector.__str__ = RVector_str
_pygimli_.CVector.__str__ = RVector_str
_pygimli_.BVector.__str__ = RVector_str
_pygimli_.IVector.__str__ = RVector_str
_pygimli_.RVector3.__str__ = RVector3_str
_pygimli_.R3Vector.__str__ = R3Vector_str

_pygimli_.RMatrix.__str__ = RMatrix_str
_pygimli_.CMatrix.__str__ = CMatrix_str
_pygimli_.Line.__str__ = Line_str
_pygimli_.Mesh.__str__ = Mesh_str
_pygimli_.DataContainer.__str__ = Data_str
#_pygimli_.stdVectorIndex.size = _pygimli_.stdVectorIndex.__len__
#_pygimli_.stdVectorIndex.__str__ = RVector_str

############################
# compatibility stuff
############################

def nonzero_test(self):
    raise BaseException("Warning! there is no 'and' and 'or' for "
                        "BVector and RVector. " +
                        "Use binary operators '&' or '|' instead. " +
                        "If you looking for the nonzero test, use len(v) > 0")

_pygimli_.RVector.__nonzero__ = nonzero_test
_pygimli_.RVector.__bool__ = nonzero_test
_pygimli_.R3Vector.__nonzero__ = nonzero_test
_pygimli_.R3Vector.__bool__ = nonzero_test
_pygimli_.BVector.__nonzero__ = nonzero_test
_pygimli_.BVector.__bool__ = nonzero_test
_pygimli_.CVector.__nonzero__ = nonzero_test
_pygimli_.CVector.__bool__ = nonzero_test
_pygimli_.IVector.__nonzero__ = nonzero_test
_pygimli_.IVector.__bool__ = nonzero_test

def __RVectorPower(self, m):
    return pow(self, m)
_pygimli_.RVector.__pow__ = __RVectorPower

############################
# Indexing [] operator for RVector, CVector,
#                          RVector3, R3Vector, RMatrix, CMatrix
############################
def __getVal(self, idx):
    """Hell slow"""
    # print("__getVal")

    if isinstance(idx, BVector) or isinstance(idx, IVector):
        return self(idx)
    elif isinstance(idx, stdVectorSIndex) or isinstance(idx, stdVectorIndex):
        return self(idx)
    elif isinstance(idx, list) or hasattr(idx, '__iter__'):
        idxL = _pygimli_.IVector(len(idx))
        for i, ix in enumerate(idx):
            if hasattr(ix, 'dtype'):
                if ix.dtype == bool:
                    idxL[i] = int(i)
                else:
                    idxL[i] = int(ix)
            else:
                idxL[i] = int(ix)

        return self(idxL)

    elif isinstance(idx, slice):
        if idx.step is None:
            s = idx.start
            e = idx.stop
            if s is None:
                s = 0
            if e is None:
                e = len(self)
# print('#'*100, s, e)
            return self.getVal(int(s), int(e))
        else:
            ids = range(idx.start, idx.stop, idx.step)
            if len(ids):
                return self(ids)
            else:
                raise Exception("slice invalid")

    elif idx == -1:
        idx = len(self) - 1

    return self.getVal(int(idx))

def __setVal(self, idx, val):

    if isinstance(idx, slice):
        if idx.step is None:
            self.setVal(float(val), int(idx.start), int(idx.stop))
            return
        else:
            print("not yet implemented")
    elif isinstance(idx, tuple):
        #print(idx, type(idx))
        self.rowR(int(idx[0])).setVal(val, int(idx[1]))
        return

    self.setVal(val, idx)

def __getValMatrix(self, idx):

    #    print(idx, type(idx))
    if isinstance(idx, slice):
        step = idx.step
        if step is None:
            step = 1
        start = idx.start
        if start is None:
            start = 0
        stop = idx.stop
        if stop is None:
            stop = len(self)

        return [self.rowR(i) for i in range(start, stop, step)]

    elif isinstance(idx, tuple):
        #print(idx, type(idx))
        if isinstance(idx[0], slice):
            if isinstance(idx[1], int):
                tmp = self.__getitem__(idx[0])
                ret = _pygimli_.RVector(len(tmp))
                for i, t in enumerate(tmp):
                    ret[i] = t[idx[1]]
                return ret
        else:
            return self.rowR(int(idx[0])).__getitem__(idx[1])

    if idx == -1:
        idx = len(self) - 1

    return self.rowR(idx)

_pygimli_.RVector.__setitem__ = __setVal
_pygimli_.RVector.__getitem__ = __getVal  # very slow -- inline is better

_pygimli_.CVector.__setitem__ = __setVal
_pygimli_.CVector.__getitem__ = __getVal  # very slow -- inline is better

_pygimli_.BVector.__setitem__ = __setVal
_pygimli_.BVector.__getitem__ = __getVal  # very slow -- inline is better

_pygimli_.IVector.__setitem__ = __setVal
_pygimli_.IVector.__getitem__ = __getVal  # very slow -- inline is better

_pygimli_.RVector3.__setitem__ = __setVal

_pygimli_.RMatrix.__getitem__ = __getValMatrix  # very slow -- inline is better
_pygimli_.RMatrix.__setitem__ = __setVal

_pygimli_.CMatrix.__getitem__ = __getValMatrix  # very slow -- inline is better
_pygimli_.CMatrix.__setitem__ = __setVal


############################
# len(RVector), RMatrix
############################
def PGVector_len(self):
    return self.size()

_pygimli_.RVector.__len__ = PGVector_len
_pygimli_.R3Vector.__len__ = PGVector_len
_pygimli_.BVector.__len__ = PGVector_len
_pygimli_.CVector.__len__ = PGVector_len
_pygimli_.IVector.__len__ = PGVector_len

def RMatrix_len(self):
    return self.rows()
_pygimli_.RMatrix.__len__ = RMatrix_len
_pygimli_.CMatrix.__len__ = RMatrix_len


############################
# Iterator support for RVector allow to apply python build-ins
############################
class VectorIter:

    def __init__(self, vec):
        self.it = vec.beginPyIter()
        self.vec = vec

    def __iter__(self):
        return self

    # this is for python < 3
    def next(self):
        return self.it.nextForPy()

    # this is the same but for python > 3
    def __next__(self):
        return self.it.nextForPy()


def __VectorIterCall__(self):
    return VectorIter(self)
    # don't use pygimli iterators here until the reference for temporary
    # vectors are collected
    # return _pygimli_.RVectorIter(self.beginPyIter())

_pygimli_.RVector.__iter__ = __VectorIterCall__
_pygimli_.R3Vector.__iter__ = __VectorIterCall__
_pygimli_.BVector.__iter__ = __VectorIterCall__
_pygimli_.IVector.__iter__ = __VectorIterCall__
_pygimli_.CVector.__iter__ = __VectorIterCall__

class DefaultContainerIter:

    def __init__(self, vec):
        self.vec = vec
        self.length = len(vec)
        self.pos = -1

    def __iter__(self):
        return self

    def next(self):
        return self.__next__()

    # this is the same but for python > 3
    def __next__(self):
        self.pos += 1
        if self.pos == self.length:
            raise StopIteration()
        else:
            return self.vec[self.pos]

def __MatIterCall__(self):
    return DefaultContainerIter(self)

_pygimli_.RMatrix.__iter__ = __MatIterCall__
_pygimli_.CMatrix.__iter__ = __MatIterCall__

class Vector3Iter (VectorIter):
    def __init__(self, vec):
        self.vec = vec
        self.length = 3
        self.pos = -1

def __Vector3IterCall__(self):
    return Vector3Iter(self)

_pygimli_.RVector3.__iter__ = __Vector3IterCall__


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
    # probably fixed!!!
    #import numpy as np
    # we need to copy the array until we can handle increasing the reference
    # counter in self.array() else it leads to strange behaviour
    # test in testRValueConverter.py:testNumpyFromRVec()
    #return np.array(self.array())
    return self.array()

_pygimli_.RVector.__array__ = __RVectorArrayCall__
_pygimli_.R3Vector.__array__ = __RVectorArrayCall__
_pygimli_.RVector3.__array__ = __RVector3ArrayCall__

# hackish until stdVectorRVector3 will be removed
def __stdVectorRVector3ArrayCall(self, idx=None):
    if idx:
        print(self)
        print(idx)
        raise Exception("we need to fix this")
    return _pygimli_.stdVectorRVector3ToR3Vector(self).array()

_pygimli_.stdVectorRVector3.__array__ = __stdVectorRVector3ArrayCall

#_pygimli_.RVector3.__array__ = _pygimli_.RVector3.array
#del _pygimli_.RVector.__array__

############################
# non automatic exposed functions
############################

def abs(v):
    if isinstance(v, _pygimli_.CVector):
        return _pygimli_.mag(v)
    elif type(v) == _pygimli_.RMatrix:
        for i in range(len(v)):
            v[i] = _pygimli_.abs(v[i])
        return v

    return _pygimli_.fabs(v)


########################################################
# compare operators for stdVector
########################################################

def __CMP_stdVectorSIndex__(self, val):
    raise Exception("__CMP_stdVectorSIndex__ do not use")
    return 0
    print("__CMP_stdVectorSIndex__")
    ret = _pygimli_.BVector(len(self))
    for i, v in enumerate(ret):
        print(self[i] < val, int(self[i] < val))
        #v = int(self[i] < val)
    print(ret)


def __EQ_stdVectorSIndex__(self, val):
    raise Exception("__EQ_stdVectorSIndex__ do not use")

    #ret = _pygimli_.BVector(len(self))
    # for i, v in enumerate(ret):
    #    print(self[i] == val, int(self[i] == val))
    #    v = self[i] == val
    # print(ret)

_pygimli_.stdVectorSIndex.__cmp__ = __CMP_stdVectorSIndex__
_pygimli_.stdVectorSIndex.__eq__ = __EQ_stdVectorSIndex__


############################
# usefull stuff
############################
def toIVector(v):
    ret = _pygimli_.IVector(len(v), 0)
    for i, r in enumerate(v):
        ret[i] = int(r)
    return ret


# DEPRECATED for backward compatibility should be removed
def asvector(array):
    return _pygimli_.RVector(array)


############################
# ???
############################
# def RVector_ArrayInit(self):
    ##self.ndim = [1, self.size()]
    #self.ndim = 1
#_pygimli_.RVector.isArray = RVector_ArrayInit


#
# Construct RVector from numpy array , opposite to asarray(RVector)
#
# def asvector(arr):
    #r = _pygimli_.RVector(len(arr))

    # for i, v in enumerate(arr):
    #r[i] = v

    # return r

# PEP conform version str with SVN revision number
def __svnversion__(path=__path__[0]):
    p = subprocess.Popen("svnversion -n %s" % path, shell=True,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdout, _) = p.communicate()
    version = str(stdout)
    if version == "Nicht versioniertes Verzeichnis":
        return ""
    else:
        return "_rev" + version

__version__ = _pygimli_.versionStr() + __svnversion__()

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

# CR: What was this good for??
#def __getRotMatrix__(self, mat):
#   getRotMatrix__(self, mat)
#_pygimli_.RQuaternion.rotMatrix = __getRotMatrix__

# some rvector helpers


def randN(n):
    """Create RVector of length n with normally distributed random numbers."""
    r = _pygimli_.RVector(n)
    _pygimli_.randn(r)
    return r

def test(show=False, coverage=False):
    """Test all code examples in docstrings using pytest."""
    try:
        import pytest
    except ImportError:
        raise ImportError("pytest is required to run test suite. " + \
                          "Try 'sudo pip install pytest'.")

    from matplotlib import pyplot as plt
    from pygimli.utils import opt_import
    pc = opt_import("pytest_cov", "create a coverage report")

    old_backend = plt.get_backend()
    if not show:
        plt.switch_backend("Agg")
    cwd = __path__[0]
    cfg = os.path.join(cwd, '../tests/setup.cfg')
    cmd = ""
    if os.path.exists(cfg):
        cmd += "-c %s " % cfg
    if pc and coverage:
        cmd += "--cov pygimli --cov-report coveralls --cov-report html " + \
               "--cov-config %s " % cfg.replace("setup.cfg", ".coveragerc")
    cmd += "%s" % cwd
    try:
        pytest.main(cmd)
    finally:
        plt.switch_backend(old_backend)
