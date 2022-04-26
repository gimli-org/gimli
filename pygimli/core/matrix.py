# -*- coding: utf-8 -*-
"""Some matrix specialization."""

import time
import numpy as np

from .core import pgcore

from .core import (CMatrix, CSparseMapMatrix, CSparseMatrix,
                   RSparseMapMatrix, RSparseMatrix, ElementMatrix,
                   IVector, MatrixBase, R3Vector, RVector)

from .logger import critical, warn, error
from .base import isArray


# make core matrices (now in pgcore, later pg.core) available here for brevity
# useful aliases
IdentityMatrix = pgcore.IdentityMatrix

Matrix = pgcore.RMatrix
SparseMatrix = pgcore.RSparseMatrix
SparseMapMatrix = pgcore.RSparseMapMatrix
BlockMatrix = pgcore.RBlockMatrix

# General Monkeypatch core classes
__Matrices = [pgcore.MatrixBase,
              #pgcore.RSparseMatrix,
              #pgcore.RSparseMapMatrix,
              #pgcore.CSparseMatrix,
              #pgcore.CSparseMapMatrix,
              #pgcore.RBlockMatrix,
              #pgcore.IdentityMatrix
              ]

for m in __Matrices:
    m.ndim = 2
    m.__len__ = lambda self: self.rows()
    m.shape = property(lambda self: (self.rows(), self.cols()))

pgcore.RMatrix.dtype = np.float
pgcore.CMatrix.dtype = np.complex
pgcore.RSparseMapMatrix.dtype = np.float
pgcore.CSparseMapMatrix.dtype = np.complex
pgcore.RSparseMatrix.dtype = np.float
pgcore.CSparseMatrix.dtype = np.complex
pgcore.RDenseMatrix.dtype = np.float


def __RMatrix_str(self):
    import pygimli as pg
    s = "RMatrix: " + str(self.rows()) + " x " + str(self.cols())

    if self.rows() < 6:
        s += '\n'

        for row in self:
            for v in row:
                s += pg.pf(v).rjust(8)
            s+='\n'
    return s

def __CMatrix_str(self):
    s = "CMatrix: " + str(self.rows()) + " x " + str(self.cols())

    if self.rows() < 6:
        s += '\n'
        for v in range(self.rows()):
            s += self[v].__str__() + '\n'
    return s


def __ElementMatrix_str(self):
    """Show entries of an ElementMatrix."""
    import pygimli as pg
    self.integrate()

    if self.mat_RM().cols() == 0 and self.mat_RM().rows() == 0:
        return 'Empty ElementMatrix\n'

    maxRowID = int(np.log10(max(self.rowIDs())))+2

    if self.multR is not None:
        if (pg.isScalar(self.multR) and self.multR != 1.0) or \
            pg.isPos(self.multR) or \
            pg.isArray(self.multR) or\
            pg.isMatrix(self.multR):
            s = '\n ' + f' multR = {self.multR} (applied)' + '\n ' + ' ' * maxRowID
        elif pg.isScalar(self.multR, 1.0):
            s = '\n ' + ' ' * maxRowID
        else:
            s = '\n ' + f' multR = {self.multR} (unknown how to apply)' + '\n ' + ' ' * maxRowID
    else:
        s = '\n ' + ' ' * maxRowID

    # print(self.mat())
    # print(self.colIDs())
    # print(self.rowIDs())
    for i in range(self.mat_RM().cols()):
        s += str(self.colIDs()[i]).rjust(9)
    s += '\n'

    s += '  ' + '-' * self.mat_RM().cols()*(9 + maxRowID) + '-\n'

    for i in range(self.mat_RM().rows()):
        s += str(self.rowIDs()[i]).rjust(maxRowID) + " :"
        if isinstance(self.multR, (int, float)):
            for v in self.row_RM(i)*self.multR:
                s += pg.pf(v).rjust(9)
        elif pg.isPos(self.multR):
            if self.mat_RM().cols() == len(self.multR):
                for v in self.row_RM(i)*self.multR:
                    s += pg.pf(v).rjust(9)
            else:
                print(self.row_RM)
                print(self.multR)
                pg.critical('invalid element multR.')
        elif pg.isArray(self.multR, self.mat_RM().cols()):
            for v in self.row_RM(i)*self.multR:
                s += pg.pf(v).rjust(9)

        elif pg.isMatrix(self.multR):
            if pg.isArray(self.multR.flatten(), self.mat_RM().cols()):
                for v in self.row_RM(i)*self.multR.flatten():
                    s += pg.pf(v).rjust(9)
            else:
                print(self.row_RM)
                print(self.multR)
                pg.critical('invalid matrix element multR.')

        elif self.multR is not None:
            print(self.mat_RM())
            print(self.multR)
            pg.critical('invalid element multR.')
        else:
            for v in self.row_RM(i):
                s += pg.pf(v).rjust(9)
        s += '\n'
    return s


pgcore.RMatrix.__repr__ = __RMatrix_str
pgcore.RDenseMatrix.__repr__ = __RMatrix_str
pgcore.CMatrix.__repr__ = __CMatrix_str
pgcore.ElementMatrix.__repr__ = __ElementMatrix_str


def __CopyRMatrixTranspose__(self):
    return pgcore.RMatrix(np.array(self).T)

pgcore.RMatrix.T = property(__CopyRMatrixTranspose__)


def __SparseMatrix_str(self):
    """Show entries of an ElementMatrix."""
    import pygimli as pg

    S = pg.utils.toSparseMapMatrix(self)
    if S.cols() == 0 and S.rows() == 0:
        return 'Empty ElementMatrix\n'

    s = "{0} size = {1} x {2}, nVals = {3}".format(type(self),
                                       S.rows(), S.cols(), S.nVals())

    if S.cols() < 20:
        s += '\n'

        M = pg.utils.toDense(self)
        for mi in M:
            for v in mi:
                if (abs(v)< 1e-12 and abs(v) > 0):
                    s += ('+'+pg.pf(v)).rjust(9)
                elif (abs(v)< 1e-12 and abs(v) < 0):
                    s += ('-'+pg.pf(v)).rjust(9)
                else:
                    s += pg.pf(v).rjust(9)

            s += '\n'
    return s

pgcore.RSparseMatrix.__repr__ =__SparseMatrix_str
pgcore.RSparseMapMatrix.__repr__ =__SparseMatrix_str


def __stdVectorRSparseMapMatrix_Mult__(self, b):
    ret = None
    if isinstance(b, list):
        # ret = np.array([pgcore.mult(self, bi) for bi in b]).T
        # print(ret)

        br = b[0]
        for i in range(1, len(b)):
            br = pgcore.cat(br, b[i])

        ret = pgcore.stdVectorR3Vector()
        pgcore.mult(self, br, ret)

    elif isinstance(b, np.ndarray) and b.ndim == 2:
        ret = pgcore.stdVectorR3Vector()
        pgcore.mult(self, b.flatten('F'), ret)

    elif isArray(b):
        ret = pgcore.stdVectorRVector()
        pgcore.mult(self, b, ret)

    else:
        print(b)
        pg.critical("Don't know how to convert b")

    return ret
pgcore.stdVectorRSparseMapMatrix.__mul__ = __stdVectorRSparseMapMatrix_Mult__


def __RVector_format(self, f):
    print(f)
    return str(self)

pgcore.RVector.__format__ = __RVector_format


# Special Monkeypatch core classes
# keep original method
__BlockMatrix_addMatrix__ = pgcore.RBlockMatrix.addMatrix

def __BlockMatrix_addMatrix_happy_GC__(self, M, row=None, col=None,
                                       scale=1.0, transpose=False):
    """Add an existing matrix to this block matrix and return a unique index.

    As long row and col are None, the Matrix will not be used until a matrix
    entry is has been added.

    Monkeypatched version to increase the reference counter of M to keep the
    garbage collector happy.

    TODO
    ----
        * Add numpy matrices or convertable
        * Transpose is only for 1d arrays. Needed for matrices?

    Parameters
    ----------
    M: pg.matrix.Matrix | pg.Vector | 1d iterable
        Matrix to add to the block.
    row: long
        Starting row index.
    col: long
        Starting column index.
    scale: float[1.0]
        Scale all matrix entries.
    transpose: bool [False]
        Transpose the matrix.
    """
    if M.ndim == 1:
        if transpose is False:
            _M = SparseMapMatrix(list(range(len(M))), [0]*len(M), M)
        else:
            warn('BlockMatrix add (transpose==True) ... Move me to core')
            _M = SparseMapMatrix([0]*len(M), list(range(len(M))), M)
        M = _M
    else:
        if transpose is True:
            if isinstance(M, pgcore.RSparseMapMatrix):
                warn('BlockMatrix add (transpose==True) ... Move me to core')
                v = pgcore.RVector()
                i = pgcore.IndexArray([0])
                j = pgcore.IndexArray([0])
                M.fillArrays(v, i, j)
                M = SparseMapMatrix(j, i, v)
            else:
                critical("don't know yet how to add transpose matrix of type",
                         type(M))

    if not hasattr(self, '__mats__'):
        self.__mats__ = []
    self.__mats__.append(M)

    matrixID = __BlockMatrix_addMatrix__(self, M)

    if row is not None and col is not None:
        self.addMatrixEntry(matrixID, row, col, scale)

    return matrixID


def __BlockMatrix_str__(self):
    string = ("pg.matrix.BlockMatrix of size %d x %d consisting of %d "
              "submatrices.")
    return string % (self.rows(), self.cols(), len(self.entries()))


pgcore.RBlockMatrix.addMatrix = __BlockMatrix_addMatrix_happy_GC__
pgcore.RBlockMatrix.add = __BlockMatrix_addMatrix_happy_GC__
pgcore.RBlockMatrix.__repr__ = __BlockMatrix_str__
pgcore.RBlockMatrix.ndim = 2
# pgcore.CBlockMatrix.addMatrix = __BlockMatrix_addMatrix_happy_GC__
# pgcore.CBlockMatrix.add = __BlockMatrix_addMatrix_happy_GC__


def __SparseMatrixEqual__(self, T):
    """Compare two SparseMatrices"""
    from pygimli.utils import sparseMatrix2Array

    if self.shape[0] != T.shape[0] or self.shape[1] != T.shape[1]:
        warn(f'Compare sizes invalid {self.shape} vs. {T.shape}: ')
        return False

    rowsA, colsA, valsA = sparseMatrix2Array(self, indices=True, getInCRS=True)
    rowsB, colsB, valsB = sparseMatrix2Array(T, indices=True, getInCRS=True)

    if len(valsA) != len(valsB):
        warn(f'Compare value sizes invalid: {len(valsA)}, {len(valsB)}')
        return False

    meanA = np.mean(abs(valsA))
    nAB = np.linalg.norm(valsA-valsB)

    if np.isnan(nAB):
        print(valsA)
        print(valsB)
        error('norm(A-B) is nan')

    if pgcore.deepDebug() == -1:
        print(self, T)
        print('diff rows:', np.linalg.norm(np.array(rowsA)-np.array(rowsB)))
        print('diff cols:', np.linalg.norm(np.array(colsA)-np.array(colsB)))

        print(f'|A-B|={nAB}', f'mean(A)={meanA} mean(B)={np.mean(abs(valsB))}')

        print(f'nAB/meanA = {nAB/meanA}')

    if meanA > 1e-10:
        return np.all(rowsA == rowsB) and np.all(colsA == colsB) and nAB/meanA < 1e-12
    else:
        return nAB < 1e-12

pgcore.RSparseMatrix.__eq__ = __SparseMatrixEqual__
pgcore.RSparseMapMatrix.__eq__ = __SparseMatrixEqual__


def __SparseMatrixCopy__(self):
    """Create a copy."""
    return pgcore.RSparseMatrix(self)

pgcore.RSparseMatrix.copy = __SparseMatrixCopy__


def __SparseMapMatrixCopy__(self):
    """Create a copy."""
    return pgcore.RSparseMapMatrix(self)

pgcore.RSparseMapMatrix.copy = __SparseMapMatrixCopy__
