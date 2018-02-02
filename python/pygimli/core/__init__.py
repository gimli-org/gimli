# -*- coding: utf-8 -*-

"""
"""

import sys
import os
import subprocess
import traceback

if sys.platform == 'win32':
    os.environ['PATH'] = __path__[0] + ';' + os.environ['PATH']

_pygimli_ = None


###############################################################################
# TEMP: Avoid import errors for py bindings built before core submodule
def bindingpath(relpath):
    path = os.path.join(os.path.dirname(__file__), relpath)
    return os.path.abspath(os.path.join(path, "_pygimli_.so"))

new = bindingpath(".")
old = bindingpath("..")

if not os.path.isfile(new) and os.path.isfile(old):
    print("INFO: Moving %s to %s" % (old, new))
    os.rename(old, new)
###############################################################################
try:
    from . import _pygimli_
    from . _pygimli_ import *
except ImportError as e:
    print(e)
    traceback.print_exc(file=sys.stdout)
    sys.stderr.write("ERROR: cannot import the library '_pygimli_'.\n")

#######################################
###  Global convinience functions #####
#######################################

#_pygimli_.load = None
from pygimli.io import load
from pygimli.viewer import show, plt, wait
from pygimli.solver import solve
from pygimli.meshtools import interpolate

def showNow():
    pass
#    showLater(0)  # not working anymore

__swatch__ = _pygimli_.Stopwatch()

def tic(msg=None):
    """Start global stopwatch."""
    if msg is not None:
        print(msg)
    __swatch__.start()


def toc(msg=None):
    """Print elapsed time since stopwatch was started."""
    if msg is not None:
        print(msg)
    print('Elapsed time is:', dur(), "s")


def dur():
    """Return time since the stopwatch was started."""
    return __swatch__.duration()


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
    return("RVector3: (" + str(self.x()) + ", " +
           str(self.y()) + ", " + str(self.z()) + ")")


def R3Vector_str(self):
    return "R3Vector: n=" + str(self.size())


def RMatrix_str(self):
    s = "RMatrix: " + str(self.rows()) + " x " + str(self.cols())

    if self.rows() < 6:
        s += '\n'
        for v in range(self.rows()):
            s += self[v].__str__(True) + '\n'
    return s


def CMatrix_str(self):
    s = "CMatrix: " + str(self.rows()) + " x " + str(self.cols())

    if self.rows() < 6:
        s += '\n'
        for v in range(self.rows()):
            s += self[v].__str__(True) + '\n'
    return s


def Line_str(self):
    return "Line: " + str(self.p0()) + "  " + str(self.p1())


def Mesh_str(self):
    return ("Mesh: Nodes: " + str(self.nodeCount()) + " Cells: " +
            str(self.cellCount()) + " Boundaries: " +
            str(self.boundaryCount()))


def Data_str(self):
    return (
        "Data: Sensors: " +
        str(self.sensorCount()) + " data: " + str(self.size())
    )


def ElementMatrix_str(self):

    s = ''
    for i in range(self.size()):
        s += str(self.idx(i)) + "\t: "

        for j in range(self.size()):
            s += str(self.getVal(i, j)) + " "
        s += '\n'
    return s

def MeshEntity_str(self):
    """Give mesh entity infos."""
    s = str(type(self))
    s += '\tID: ' + str(self.id()) + \
         ', Marker: ' + str(self.marker()) + \
         ', Size: ' + str(self.size()) + '\n'
    for n in self.nodes():
        s += '\t' + str(n.id()) + " " + str(n.pos()) + "\n"
    return s


_pygimli_.RVector.__str__ = RVector_str
_pygimli_.CVector.__str__ = RVector_str
_pygimli_.BVector.__str__ = RVector_str
_pygimli_.IVector.__str__ = RVector_str
_pygimli_.IndexArray.__str__ = RVector_str
_pygimli_.RVector3.__str__ = RVector3_str
_pygimli_.R3Vector.__str__ = R3Vector_str

_pygimli_.RMatrix.__str__ = RMatrix_str
_pygimli_.CMatrix.__str__ = CMatrix_str
_pygimli_.Line.__str__ = Line_str
_pygimli_.Mesh.__str__ = Mesh_str
_pygimli_.DataContainer.__str__ = Data_str
_pygimli_.ElementMatrix.__str__ = ElementMatrix_str
_pygimli_.MeshEntity.__str__ = MeshEntity_str
# _pygimli_.stdVectorIndex.size = _pygimli_.stdVectorIndex.__len__
# _pygimli_.stdVectorIndex.__str__ = RVector_str

############################
# compatibility stuff
############################


def nonzero_test(self):
    raise BaseException("Warning! there is no 'and' and 'or' for "
                        "BVector and RVector. " +
                        "Use binary operators '&' or '|' instead. " +
                        "If you looking for the nonzero test, use len(v) > 0")

def np_round__(self, r):
    return np.round(self.array(), r)

_pygimli_.RVector.__nonzero__ = nonzero_test
_pygimli_.RVector.__bool__ = nonzero_test
_pygimli_.RVector.__round__ = np_round__
_pygimli_.R3Vector.__nonzero__ = nonzero_test
_pygimli_.R3Vector.__bool__ = nonzero_test
_pygimli_.BVector.__nonzero__ = nonzero_test
_pygimli_.BVector.__bool__ = nonzero_test
_pygimli_.CVector.__nonzero__ = nonzero_test
_pygimli_.CVector.__bool__ = nonzero_test
_pygimli_.IVector.__nonzero__ = nonzero_test
_pygimli_.IVector.__bool__ = nonzero_test
_pygimli_.IndexArray.__nonzero__ = nonzero_test
_pygimli_.IndexArray.__bool__ = nonzero_test

######################
# special constructors
######################

# Overwrite constructor for IndexArray
# This seams ugly but necessary until we can recognize numpy array in
# custom_rvalue
__origIndexArrayInit__ = _pygimli_.IndexArray.__init__


def __newIndexArrayInit__(self, arr, val=None):
    """"""
    # print("Custom IndexArray", arr, val)
    if hasattr(arr, 'dtype') and hasattr(arr, '__iter__'):
        __origIndexArrayInit__(self, [int(a) for a in arr])
    else:
        if val:
            __origIndexArrayInit__(self, arr, val)
        else:
            __origIndexArrayInit__(self, arr)

_pygimli_.IndexArray.__init__ = __newIndexArrayInit__

# Overwrite constructor for BVector
# This seams ugly but necessary until we can recognize numpy array in
# custom_rvalue
__origBVectorInit__ = _pygimli_.BVector.__init__


def __newBVectorInit__(self, arr, val=None):
    if hasattr(arr, 'dtype') and hasattr(arr, '__iter__'):
        # this is hell slow .. better in custom_rvalue.cpp or in
        # vector.h directly from pyobject
        __origBVectorInit__(self, len(arr))
        for i, a in enumerate(arr):
            self.setVal(bool(a), i)
    else:
        if val:
            __origBVectorInit__(self, arr, val)
        else:
            __origBVectorInit__(self, arr)

_pygimli_.BVector.__init__ = __newBVectorInit__


######################
# special overwrites
######################

# RVector + int fails .. so we need to tweak this command
__oldRVectorAdd__ = _pygimli_.RVector.__add__
def __newRVectorAdd__(a, b):
    if type(b) == int:
        return __oldRVectorAdd__(a, float(b))
    if type(a) == int:
        return __oldRVectorAdd__(float(a), b)
    return __oldRVectorAdd__(a, b)
_pygimli_.RVector.__add__ = __newRVectorAdd__

__oldRVectorSub__ = _pygimli_.RVector.__sub__
def __newRVectorSub__(a, b):
    if type(b) == int:
        return __oldRVectorSub__(a, float(b))
    if type(a) == int:
        return __oldRVectorSub__(float(a), b)
    return __oldRVectorSub__(a, b)
_pygimli_.RVector.__sub__ = __newRVectorSub__

__oldRVectorMul__ = _pygimli_.RVector.__mul__
def __newRVectorMul__(a, b):
    if type(b) == int:
        return __oldRVectorMul__(a, float(b))
    if type(a) == int:
        return __oldRVectorMul__(float(a), b)
    return __oldRVectorMul__(a, b)
_pygimli_.RVector.__mul__ = __newRVectorMul__

try:
    __oldRVectorTrueDiv__ = _pygimli_.RVector.__truediv__
    def __newRVectorTrueDiv__(a, b):
        if type(b) == int:
            return __oldRVectorTrueDiv__(a, float(b))
        if type(a) == int:
            return __oldRVectorTrueDiv__(float(a), b)
        return __oldRVectorTrueDiv__(a, b)
    _pygimli_.RVector.__truediv__ = __newRVectorTrueDiv__
except:
    __oldRVectorTrueDiv__ = _pygimli_.RVector.__div__
    def __newRVectorTrueDiv__(a, b):
        if type(b) == int:
            return __oldRVectorTrueDiv__(a, float(b))
        if type(a) == int:
            return __oldRVectorTrueDiv__(float(a), b)
        return __oldRVectorTrueDiv__(a, b)
    _pygimli_.RVector.__div__ = __newRVectorTrueDiv__


################################################################################
# override wrong default conversion from int to IndexArray(int) for setVal     #
################################################################################
__origRVectorSetVal__ = _pygimli_.RVector.setVal
def __newRVectorSetVal__(self, *args, **kwargs):
    #print('__newRVectorSetVal__', *args, **kwargs)
    if len(args)==2:
        if isinstance(args[1], int):
            if args[1] < 0:
                return __origRVectorSetVal__(self, args[0], i=len(self)+args[1])
            else:
                return __origRVectorSetVal__(self, args[0], i=args[1])
        if isinstance(args[1], _pygimli_.BVector):
            return __origRVectorSetVal__(self, args[0], bv=args[1])
    return __origRVectorSetVal__(self, *args, **kwargs)
_pygimli_.RVector.setVal = __newRVectorSetVal__


__origR3VectorSetVal__ = _pygimli_.R3Vector.setVal
def __newR3VectorSetVal__(self, *args, **kwargs):
    #print('__newRVectorSetVal__', *args, **kwargs)
    if len(args)==2:
        if isinstance(args[1], int):
            return __origR3VectorSetVal__(self, args[0], i=args[1])
        if isinstance(args[1], _pygimli_.BVector):
            return __origR3VectorSetVal__(self, args[0], bv=args[1])
    return __origR3VectorSetVal__(self, *args, **kwargs)
_pygimli_.R3Vector.setVal = __newR3VectorSetVal__


__origBVectorSetVal__ = _pygimli_.BVector.setVal
def __newBVectorSetVal__(self, *args, **kwargs):
    if len(args)==2:
        if isinstance(args[1], int):
            return __origBVectorSetVal__(self, args[0], i=args[1])
        if isinstance(args[1], _pygimli_.BVector):
            return __origBVectorSetVal__(self, args[0], bv=args[1])
    return __origBVectorSetVal__(self, *args, **kwargs)
_pygimli_.BVector.setVal = __newBVectorSetVal__

__origCVectorSetVal__ = _pygimli_.CVector.setVal
def __newCVectorSetVal__(self, *args, **kwargs):
    if len(args)==2:
        if isinstance(args[1], int):
            return __origCVectorSetVal__(self, args[0], i=args[1])
        if isinstance(args[1], _pygimli_.BVector):
            return __origCVectorSetVal__(self, args[0], bv=args[1])
    return __origCVectorSetVal__(self, *args, **kwargs)
_pygimli_.CVector.setVal = __newCVectorSetVal__

__origIVectorSetVal__ = _pygimli_.IVector.setVal
def __newIVectorSetVal__(self, *args, **kwargs):
    if len(args)==2:
        if isinstance(args[1], int):
            return __origIVectorSetVal__(self, args[0], i=args[1])
        if isinstance(args[1], _pygimli_.BVector):
            return __origIVectorSetVal__(self, args[0], bv=args[1])
    return __origIVectorSetVal__(self, *args, **kwargs)
_pygimli_.IVector.setVal = __newIVectorSetVal__

__origIndexArraySetVal__ = _pygimli_.IndexArray.setVal
def __newIndexArraySetVal__(self, *args, **kwargs):
    if len(args)==2:
        if isinstance(args[1], int):
            return __origIndexArraySetVal__(self, args[0], i=args[1])
        if isinstance(args[1], _pygimli_.BVector):
            return __origIndexArraySetVal__(self, args[0], bv=args[1])
    return __origIndexArraySetVal__(self, *args, **kwargs)
_pygimli_.IndexArray.setVal = __newIndexArraySetVal__


############################
# Indexing [] operator for RVector, CVector,
#                          RVector3, R3Vector, RMatrix, CMatrix
############################
def __getVal(self, idx):
    """Hell slow"""

    # print("getval", type(idx), idx)
    # print(dir(idx))
    if isinstance(idx, _pygimli_.BVector) or isinstance(idx, _pygimli_.IVector) or isinstance(idx, _pygimli_.IndexArray):
        # print("BVector, IVector, IndexArray", idx)
        return self(idx)
    elif isinstance(idx, _pygimli_.stdVectorSIndex): # // or isinstance(idx, stdVectorIndex):
        # print("stdVectorSIndex", idx)
        return self(idx)
    elif isinstance(idx, slice):

        s = idx.start
        e = idx.stop
        if s is None:
            s = 0
        if e is None:
            e = len(self)

        if idx.step is None:
            return self.getVal(int(s), int(e))
        else:
            #print(s,e,idx.step)
            step = idx.step
            if step < 0 and idx.start is None and idx.stop is None:
                ids = range(e-1, s-1, idx.step)
            else:
                ids = range(s, e, idx.step)

            if len(ids):
                return self(ids)
            else:
                return self(0)
                #raise Exception("slice invalid")

    elif isinstance(idx, list) or hasattr(idx, '__iter__'):
        if isinstance(idx[0], int):
            return self(idx)
        elif hasattr(idx[0], 'dtype'):
            # print("numpy: ", idx[0].dtype.str, idx[0].dtype ,type(idx[0]))
            if idx[0].dtype == 'bool':
                return self([i for i, x in enumerate(idx) if x])
                # return self[np.nonzero(idx)[0]]
        # print("default")
        return self([int(a) for a in idx])

    elif idx == -1:
        idx = len(self) - 1

    return self.getVal(int(idx))


def __setVal(self, idx, val):
    """
        Index write access []
    """
    # print("__setVal", self, 'idx', idx, 'val:', val)
    if isinstance(idx, slice):
        if idx.step is None:
            self.setVal(val, int(idx.start), int(idx.stop))
            return
        else:
            print("not yet implemented")
    elif isinstance(idx, tuple):
        # print(idx, type(idx))
        self.rowR(int(idx[0])).setVal(val, int(idx[1]))
        return
    # if isinstance(idx, _pygimli_.BVector):
        # print("__setVal", self, idx, 'val:', val)
        # self.setVal(val, bv=idx)
        # return
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
        # print(idx, type(idx))
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

_pygimli_.R3Vector.__setitem__ = __setVal
_pygimli_.R3Vector.__getitem__ = __getVal  # very slow -- inline is better

_pygimli_.IndexArray.__setitem__ = __setVal
_pygimli_.IndexArray.__getitem__ = __getVal  # very slow -- inline is better

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
_pygimli_.IndexArray.__len__ = PGVector_len


def RMatrix_len(self):
    return self.rows()
_pygimli_.RMatrix.__len__ = RMatrix_len
_pygimli_.CMatrix.__len__ = RMatrix_len


_pygimli_.RVector.ndim = 1
_pygimli_.BVector.ndim = 1
_pygimli_.CVector.ndim = 1
_pygimli_.RVector3.ndim = 1
_pygimli_.IVector.ndim = 1
_pygimli_.IndexArray.ndim = 1
_pygimli_.RMatrix.ndim = 2
_pygimli_.R3Vector.ndim = 2
_pygimli_.stdVectorRVector3 = 2

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
_pygimli_.IndexArray.__iter__ = __VectorIterCall__
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


class Vector3Iter():
    """
    Simple iterator for RVector3, cause it lacks the core function
    .beginPyIter()
    """
    def __init__(self, vec):
        self.vec = vec
        self.length = 3
        self.pos = -1

    def __iter__(self):
        return self

    def next(self):
        return self.__next__()

    def __next__(self):
        self.pos += 1
        if self.pos == self.length:
            raise StopIteration()
        else:
            return self.vec[self.pos]


def __Vector3IterCall__(self):
    return Vector3Iter(self)

_pygimli_.RVector3.__iter__ = __Vector3IterCall__


# ######### c to python converter ######
# default converter from RVector3 to numpy array
def __RVector3ArrayCall__(self, dtype=None):
    #if idx:
        #print(self)
        #print(idx)
        #raise Exception("we need to fix this")
    import numpy as np
    return np.array([self.getVal(0), self.getVal(1), self.getVal(2)])

# default converter from RVector to numpy array


def __RVectorArrayCall__(self, dtype=None):

    #if idx and not isinstance(idx, numpy.dtype):
        #print("self:", self)
        #print("idx:", idx, type(idx) )
        #raise Exception("we need to fix this")
    # probably fixed!!!
    # import numpy as np
    # we need to copy the array until we can handle increasing the reference
    # counter in self.array() else it leads to strange behaviour
    # test in testRValueConverter.py:testNumpyFromRVec()
    # return np.array(self.array())
    return self.array()

# default converter from RVector to numpy array

_pygimli_.RVector.__array__ = __RVectorArrayCall__
# not yet ready handmade_wrappers.py
_pygimli_.BVector.__array__ = __RVectorArrayCall__
# not yet ready handmade_wrappers.py
# _pygimli_.IndexArray.__array__ = __RVectorArrayCall__
_pygimli_.R3Vector.__array__ = __RVectorArrayCall__
_pygimli_.RVector3.__array__ = __RVector3ArrayCall__

# hackish until stdVectorRVector3 will be removed
def __stdVectorRVector3ArrayCall(self, dtype=None):
    #if idx is not None:
        #print(self)
        #print(idx)
        #raise Exception("we need to fix this")
    return _pygimli_.stdVectorRVector3ToR3Vector(self).array()

_pygimli_.stdVectorRVector3.__array__ = __stdVectorRVector3ArrayCall

# _pygimli_.RVector3.__array__ = _pygimli_.RVector3.array
# del _pygimli_.RVector.__array__



##################################
# custom rvalues for special cases
##################################

def find(v):
    if hasattr(v, 'dtype') and hasattr(v, '__iter__'):
        # print('new find', v, _pygimli_.BVector(v))
        return _pygimli_.find(_pygimli_.BVector(v))
    else:
        # print('orig find')
        return _pygimli_.find(v)

def pow(v, p):
    """
        pow(v, int) is misinterpreted as pow(v, rvec(int))
        so we need to fix this
    """
    if isinstance(p, int):
        return _pygimli_.pow(v, float(p))
    return _pygimli_.pow(v, p)


def __RVectorPower(self, m):
    return pow(self, m)
_pygimli_.RVector.__pow__ = __RVectorPower


##################################
# usefull aliases
##################################

BlockMatrix = _pygimli_.RBlockMatrix
SparseMapMatrix = _pygimli_.RSparseMapMatrix
SparseMatrix = _pygimli_.RSparseMatrix
Vector = _pygimli_.RVector
Matrix = _pygimli_.RMatrix
Inversion = _pygimli_.RInversion
Pos = _pygimli_.RVector3

Trans = _pygimli_.RTrans
TransLinear = _pygimli_.RTransLinear
TransLin = _pygimli_.RTransLin
TransPower =_pygimli_.RTransPower
TransLog = _pygimli_.RTransLog
TransLogLU = _pygimli_.RTransLogLU
TransCotLU = _pygimli_.RTransCotLU
TransCumulative = _pygimli_.RTransCumulative


############################
# non automatic exposed functions
############################


def abs(v):
    if isinstance(v, _pygimli_.CVector):
        return _pygimli_.mag(v)
    elif isinstance(v, _pygimli_.R3Vector):
        return _pygimli_.absR3(v)
    elif isinstance(v, np.ndarray):
        return _pygimli_.absR3(v)
    elif type(v) == _pygimli_.RMatrix:
        raise BaseException("IMPLEMENTME")
        for i in range(len(v)):
            v[i] = _pygimli_.abs(v[i])
        return v

    return _pygimli_.fabs(v)


########################################################
# compare operators for stdVector
########################################################

def __CMP_stdVectorSIndex__(self, val):
    raise BaseException("__CMP_stdVectorSIndex__ do not use")

    # print("__CMP_stdVectorSIndex__")
    # ret = _pygimli_.BVector(len(self))
    # for i, v in enumerate(ret):
    #     print(self[i] < val, int(self[i] < val))
    #     #v = int(self[i] < val)
    # print(ret)


def __EQ_stdVectorSIndex__(self, val):
    raise Exception("__EQ_stdVectorSIndex__ do not use")

    # ret = _pygimli_.BVector(len(self))
    # for i, v in enumerate(ret):
    #    print(self[i] == val, int(self[i] == val))
    #    v = self[i] == val
    # print(ret)

_pygimli_.stdVectorSIndex.__cmp__ = __CMP_stdVectorSIndex__
_pygimli_.stdVectorSIndex.__eq__ = __EQ_stdVectorSIndex__

# default BVector operator == (RVector, int) will be casted to
# BVector operator == (RVector, RVector(int)) and fails
# this needs a monkey patch for BVector operator == (RVector, int)
_pygimli_.__EQ_RVector__ = _pygimli_.RVector.__eq__

def __EQ_RVector__(self, val):
    if type(val) == int:
        val = float(val)
    return _pygimli_.__EQ_RVector__(self, val)

_pygimli_.RVector.__eq__ = __EQ_RVector__

############################
# usefull stuff
############################
def toIVector(v):
    print("do not use toIVector(v) use ndarray directly .. "
          "this method will be removed soon")
    ret = _pygimli_.IVector(len(v), 0)
    for i, r in enumerate(v):
        ret[i] = int(r)
    return ret


#__catOrig__ = _pygimli_.cat

#def __cat__(v1, v2):
    #print("mycat")
    #if isinstance(v1, ndarray) and isinstance(v2, ndarray):
        #return pg.cat(RVector(v1), v2)
    #else:
        #return __catOrig__(v1, v2)

#_pygimli_.cat = __cat__



# DEPRECATED for backward compatibility should be removed
def asvector(array):
    print("do not use asvector(ndarray) use ndarray directly .. "
          "this method will be removed soon")
    return _pygimli_.RVector(array)


def __gitversion__(path=__path__[0]):
    """Return version str based on git tags and commits."""
    try:
        v = subprocess.check_output(['git', '-C', path, 'describe', '--always',
                                     '--tags', '--dirty=-modified', '--long'],
                                    stderr=subprocess.STDOUT)
        return v.decode('ascii').strip()
    except BaseException as e:
        return "unknown", e

# __version__ = __gitversion__()
# if __version__ == "unknown":
#     __version__ = _pygimli_.versionStr()


# ##########################
# We want ModellingBase with multi threading jacobian brute force
# ##########################

def __GLOBAL__response_mt_shm_(fop, model, shm, i):
    resp = fop.response_mt(model, i)

    for j in range(len(resp)):
        shm[j] = resp[j]


def __ModellingBase__createJacobian_mt__(self, model, resp):
    from math import ceil
    from multiprocessing import Process, Array
    import numpy as np

    nModel = len(model)
    nData = len(resp)

    fak = 1.05

    dModel = _pygimli_.RVector(len(model))
    nProcs = self.multiThreadJacobian()

    shm = []

    oldBertThread = self.threadCount()
    self.setThreadCount(1)

    # print("Model/Data/nProcs", nModel, nData, nProcs, int(ceil(float(nModel)/nProcs)))
    for pCount in range(int(ceil(float(nModel)/nProcs))):
        procs = []
        #if self.verbose():
        tic()
        if self.verbose():
            print("Jacobian MT:(", pCount*nProcs, "--" , (pCount+1)*nProcs, ") /" ,
              nModel, '... ')

        for i in range(int(pCount*nProcs), int((pCount+1)*nProcs)):
            if i < nModel:
                modelChange = _pygimli_.RVector(model)
                modelChange[i] *= fak
                dModel[i] = modelChange[i]-model[i]

                shm.append(Array('d', len(resp)))
                procs.append(Process(target=__GLOBAL__response_mt_shm_,
                                     args=(self, modelChange, shm[i], i)))

        for i, p in enumerate(procs):
            p.start()

        for i, p in enumerate(procs):
            p.join()

        if self.verbose():
            print(dur(), 's')
    self.setThreadCount(oldBertThread)

    for i in range(nModel):
        dData = np.array(shm[i]) - resp
        self._J.setCol(i, dData/dModel[i])


def __ModellingBase__responses_mt__(self, models, respos):

    nModel = len(models)
    nProcs = self.multiThreadJacobian()

    if nProcs == 1:
        for i, m in enumerate(models):
            respos[i] = self.response_mt(m, i)
        return

    from math import ceil
    from multiprocessing import Process, Array
    import numpy as np

    if models.ndim != 2:
        raise BaseException("models need to be a matrix(N, nModel):" +
                            str(models.shape))
    if respos.ndim != 2:
        raise BaseException("respos need to be a matrix(N, nData):" +
                            str(respos.shape))

    nData = len(respos[0])
    shm = []

    oldBertThread = self.threadCount()
    self.setThreadCount(1)

    # print("*"*100)
    # print(nModel, nProcs)
    # print("*"*100)
    for pCount in range(int(ceil(nModel/nProcs))):
        procs = []
        if self.verbose():
            print(pCount*nProcs, "/", nModel)
        for i in range(int(pCount*nProcs), int((pCount+1)*nProcs)):

            if i < nModel:
                shm.append(Array('d', nData))
                procs.append(Process(target=__GLOBAL__response_mt_shm_,
                                     args=(self, models[i], shm[i], i)))

        for i, p in enumerate(procs):
            p.start()

        for i, p in enumerate(procs):
            p.join()

    self.setThreadCount(oldBertThread)

    for i in range(nModel):
        resp = np.array(shm[i])
        respos[i] = resp


class ModellingBaseMT__(_pygimli_.ModellingBase):
    def __init__(self, mesh=None, dataContainer=None, verbose=False):
        if mesh and dataContainer:
            _pygimli_.ModellingBase.__init__(self, mesh=mesh,
                                             dataContainer=dataContainer,
                                             verbose=verbose)
        elif type(mesh) is _pygimli_.Mesh:
            _pygimli_.ModellingBase.__init__(self, mesh=mesh, verbose=verbose)
        elif dataContainer:
            _pygimli_.ModellingBase.__init__(self, dataContainer=dataContainer,
                                             verbose=verbose)
        else:
            _pygimli_.ModellingBase.__init__(self, verbose=verbose)

        self._J = _pygimli_.RMatrix()
        self.setJacobian(self._J)

ModellingBaseMT__.createJacobian_mt = __ModellingBase__createJacobian_mt__
ModellingBaseMT__.responses = __ModellingBase__responses_mt__

ModellingBase = ModellingBaseMT__


###########################
# unsorted stuff
###########################

# DEPRECATED
# _pygimli_.interpolate = _pygimli_.interpolate_GILsave__

from .matrix import *

############################
# some backward compatibility
############################

def __MeshGetCellMarker__(self):
    deprecated(msg='Mesh::cellMarker()', hint='Mesh::cellMarkers()')
    return self.cellMarkers()


def __MeshSetCellMarker__(self, m):
    deprecated(msg='Mesh::setCellMarker()', hint='Mesh::setCellMarkers()')
    return self.setCellMarkers(m)


_pygimli_.Mesh.cellMarker = __MeshGetCellMarker__
_pygimli_.Mesh.setCellMarker = __MeshSetCellMarker__
