# -*- coding: utf-8 -*-
"""Some matrix specialization."""

import time
import numpy as np
#from pygimli.core import _pygimli_ as pg

from . import _pygimli_ as pgcore
from . import (CMatrix, CSparseMapMatrix, CSparseMatrix,
               RSparseMapMatrix, RSparseMatrix, ElementMatrix,
               IVector, MatrixBase, R3Vector, RVector)

from .logger import critical, warn

# make core matrices (now in pgcor, later pg.core) available here for brevity
## Usefull Aliases
IdentityMatrix = pgcore.IdentityMatrix

Matrix = pgcore.RMatrix
SparseMatrix = pgcore.RSparseMatrix
SparseMapMatrix = pgcore.RSparseMapMatrix
BlockMatrix = pgcore.RBlockMatrix

## General Monkeypatch core classes
__Matrices = [pgcore.MatrixBase,
              pgcore.RSparseMatrix,
              pgcore.RSparseMapMatrix,
              pgcore.CSparseMatrix,
              pgcore.CSparseMapMatrix,
              pgcore.RBlockMatrix,
              pgcore.IdentityMatrix]

def __Matrix_len(self):
    return self.rows()

@property
def __MatrixShapePropery__(self):
    return (self.rows(), self.cols())

for m in __Matrices:
    m.ndim = 2
    m.__len__ = __Matrix_len
    m.shape = __MatrixShapePropery__

pgcore.RMatrix.dtype = np.float
pgcore.CMatrix.dtype = np.complex

def __RMatrix_str(self):
    s = "RMatrix: " + str(self.rows()) + " x " + str(self.cols())

    if self.rows() < 6:
        s += '\n'
        for v in range(self.rows()):
            s += self[v].__str__(True) + '\n'
    return s


def __CMatrix_str(self):
    s = "CMatrix: " + str(self.rows()) + " x " + str(self.cols())

    if self.rows() < 6:
        s += '\n'
        for v in range(self.rows()):
            s += self[v].__str__(True) + '\n'
    return s

def __ElementMatrix_str(self):
    """Show entries of an ElementMatrix."""
    import pygimli as pg
    if self.mat().cols() == 0 and self.mat().rows() == 0:
        return 'Empty ElementMatrix\n'

    maxRowID = int(np.log10(max(self.rowIDs())))+2

    s = '\n ' + ' ' * maxRowID
    # print(self.mat())
    # print(self.colIDs())
    # print(self.rowIDs())
    for i in range(self.mat().cols()):
        s += str(self.colIDs()[i]).rjust(9)
    s += '\n'

    s += '  ' + '-'*self.mat().cols()*(9 + maxRowID) + '-\n'

    for i in range(self.mat().rows()):
        s += str(self.rowIDs()[i]).rjust(maxRowID) + " :"
        for v in self.row(i):
            s += pg.pf(v).rjust(9)
        s += '\n'
    return s

pgcore.RMatrix.__str__ = __RMatrix_str
pgcore.CMatrix.__str__ = __CMatrix_str
pgcore.ElementMatrix.__str__ = __ElementMatrix_str


## Special Monkeypatch core classes
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
    from pygimli.utils import sparseMatrix2Array

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
pgcore.RBlockMatrix.__str__ = __BlockMatrix_str__
pgcore.RBlockMatrix.ndim = 2
# pgcore.CBlockMatrix.addMatrix = __BlockMatrix_addMatrix_happy_GC__
# pgcore.CBlockMatrix.add = __BlockMatrix_addMatrix_happy_GC__


def __SparseMatrixEqual__(self, T):
    """Compare two SparseMatrices"""
    from pygimli.utils import sparseMatrix2Array
    if self.rows() != T.rows() or self.cols() != T.cols():
        warn("Compare sizes invalid {0},{1} vs. {2},{3}: ".format(
            self.rows(), self.cols(), T.rows(), T.cols()))
        return False


    rowsA, colsA, valsA = sparseMatrix2Array(self, indices=True)
    rowsB, colsB, valsB = sparseMatrix2Array(T, indices=True)

    if len(valsA) != len(valsB):
        warn("Compare value sizes invalid: ", len(valsA), len(valsB))
        return False

    # print(np.linalg.norm(valsA-valsB))

    return rowsA == rowsB and \
           colsA == colsB and \
               np.linalg.norm(valsA-valsB) < 1e-14

pgcore.RSparseMatrix.__eq__ = __SparseMatrixEqual__
pgcore.RSparseMapMatrix.__eq__ = __SparseMatrixEqual__


class MultMatrix(pgcore.MatrixBase):
    """Base Matrix class for all matrix types holding a matrix."""
    def __init__(self, A, verbose=False):
        self._A = A
        self.ndim = self._A.ndim
        super(MultMatrix, self).__init__(verbose)

    @property
    def A(self):
        return self._A

    @A.setter
    def A(self, A):
        self._A = A

    def rows(self):
        """Return number of rows (using underlying matrix)."""
        return self.A.rows()  # this should be _A

    def cols(self):
        """Return number of columns (using underlying matrix)."""
        return self.A.cols()  # this should be _A

    def save(self, filename):
        """So it can be used in inversion with dosave flag"""
        pass


class MultLeftMatrix(MultMatrix):
    """Matrix consisting of actual RMatrix and lef-side vector."""

    def __init__(self, A, left, verbose=False):
        """Constructor saving matrix and vector."""
        if A.rows() != len(left):
            raise Exception("Matrix columns do not fit vector length!")
        super(MultLeftMatrix, self).__init__(A, verbose)

        self._l = left

    @property
    def l(self):  # better use left and right instead (pylint E743)?
        return self._l

    @l.setter
    def r(self, l):
        self._l = l

    def mult(self, x):
        """Multiplication from right-hand-side (dot product A*x)."""
        return self.A.mult(x) * self.l

    def transMult(self, x):
        """Multiplication from right-hand-side (dot product A.T * x)"""
        return self.A.transMult(x * self.l)

LMultRMatrix = MultLeftMatrix  # alias for backward compatibility


class MultRightMatrix(MultMatrix):
    """Some Matrix, multiplied with a right hand side vector r."""

    def __init__(self, A, r=None, verbose=False):
        super(MultRightMatrix, self).__init__(A, verbose)

        if r is None:
            self._r = pgcore.RVector(self.cols(), 1.0)
        else:
            self._r = r

    @property
    def r(self):
        return self._r

    @r.setter
    def r(self, r):
        self._r = r

    def mult(self, x):
        """Return M*x = A*(r*x)"""
        if hasattr(x, '__len__') and hasattr(self.r, '__len__'):
            if len(x) != len(self.r):
                # assuming A was complex
                # warn('need to double x')
                # print('mult:', self.A.rows(), " x " , self.A.cols(),
                #        'x:', len(x), 'r:', len(self.r))
                # print(self.perm)
                return self.A.mult(x[self.perm] * self.r)
                # return self.A.mult(pgcore.cat(x, x) * self.r)
        return self.A.mult(x * self.r)

    def transMult(self, x):
        """Return M.T*x=(A.T*x)*r"""
        # print('transmult', self.A.rows(), " x " , self.A.cols(), x, self.r, )
        return self.A.transMult(x) * self.r

RMultRMatrix = MultRightMatrix  # alias for backward compatibility


class MultLeftRightMatrix(MultMatrix):
    """Matrix consisting of actual RMatrix and left-hand-side vector."""

    def __init__(self, A, left, right, verbose=False):
        """Constructor saving matrix and vector."""
        if A.cols() != len(right):
            raise Exception("Matrix columns do not fit right vector length!")
        if A.rows() != len(left):
            raise Exception("Matrix rows do not fit left vector length!")

        super(MultLeftRightMatrix, self).__init__(A, verbose)
        self._r = right
        self._l = left

    @property
    def l(self):
        return self._l

    @l.setter
    def l(self, l):
        self._l = l

    @property
    def r(self):
        return self._r

    @r.setter
    def r(self, r):
        self._r = r

    def mult(self, x):
        """Multiplication from right-hand-side (dot product A*x)."""
        return self.A.mult(x * self._r) * self._l

    def transMult(self, x):
        """Multiplication from right-hand-side (dot product A.T*x)."""
        return self.A.transMult(x * self._l) * self._r

LRMultRMatrix = MultLeftRightMatrix  # alias for backward compatibility


class Add2Matrix(pgcore.MatrixBase):
    """Matrix by adding two matrices."""

    def __init__(self, A, B):
        super().__init__()
        self.A = A
        self.B = B
        assert A.rows() == B.rows()
        assert A.cols() == B.cols()

    def mult(self, x):
        """Return M*x = A*(r*x)"""
        return self.A.mult(x) + self.B.mult(x)

    def transMult(self, x):
        """Return M.T*x=(A.T*x)*r"""
        return self.A.transMult(x) + self.B.transMult(x)

    def cols(self):
        """Number of columns."""
        return self.A.cols()

    def rows(self):
        """Number of rows."""
        return self.A.rows()


class Mult2Matrix(pgcore.MatrixBase):
    """Matrix by multiplying two matrices.
        M*x = A * (B*x)
        M.T*x = (A*x) * B
    """

    def __init__(self, A, B):
        super().__init__()
        self.A = A
        self.B = B
        assert A.cols() == B.rows()

    def mult(self, x):
        """Return M*x = A*(B*x)"""
        return self.A.mult(self.B.mult(x))

    def transMult(self, x):
        """Return M.T*x=(A.T*x)*B"""
        return self.B.transMult(self.A.transMult(x))

    def cols(self):
        """Number of columns."""
        return self.B.cols()

    def rows(self):
        """Number of rows."""
        return self.A.rows()


class DiagonalMatrix(pgcore.MatrixBase):
    """Square matrix with a vector on the main diagonal."""

    def __init__(self, d):
        super().__init__()
        self.d = d

    def mult(self, x):
        """Return M*x = d*x (element-wise)"""
        return x * self.d

    def transMult(self, x):
        """Return M.T*x = M*x"""
        return x * self.d

    def cols(self):
        """Number of columns (length of diagonal)."""
        return len(self.d)

    def rows(self):
        """Number of rows (length of diagonal)."""
        return len(self.d)


class Cm05Matrix(pgcore.MatrixBase):
    """Matrix implicitly representing the inverse square-root."""

    def __init__(self, A, verbose=False):
        """Constructor saving matrix and vector.

        Parameters
        ----------
        A : ndarray
            numpy type (full) matrix
        """
        from scipy.linalg import eigh  # , get_blas_funcs

        if A.shape[0] != A.shape[1]:  # rows/cols for pgcore matrix
            raise Exception("Matrix must by square (and symmetric)!")

        self.size = A.shape[0]
        t = time.time()
        self.ew, self.EV = eigh(A)
        self.mul = np.sqrt(1./self.ew)
        if verbose:
            pgcore.info('(C) Time for eigenvalue decomposition:{:.1f}s'.format(
                time.time() - t))

        self.A = A
        super().__init__(verbose)  # only in Python 3

    def rows(self):
        """Return number of rows (using underlying matrix)."""
        return self.size

    def cols(self):
        """Return number of columns (using underlying matrix)."""
        return self.size

    def mult(self, x):
        """Multiplication from right-hand side (dot product)."""
        part1 = (np.dot(np.transpose(x), self.EV).T*self.mul).reshape(-1, 1)
        return self.EV.dot(part1).reshape(-1,)
#        return self.EV.dot((x.T.dot(self.EV)*self.mul).T)

    def transMult(self, x):
        """Multiplication from right-hand side (dot product)."""
        return self.mult(x)  # matrix is symmetric by definition


class NDMatrix(BlockMatrix):
    """Diagonal block (block-Jacobi) matrix derived from pg.matrix.BlockMatrix.

    (to be moved to a better place at a later stage)
    """

    def __init__(self, num, nrows, ncols):
        super(NDMatrix, self).__init__()  # call inherited init function
        self.Ji = []  # list of individual block matrices
        for i in range(num):
            self.Ji.append(pgcore.Matrix())
            self.Ji[-1].resize(nrows, ncols)
            n = self.addMatrix(self.Ji[-1])
            self.addMatrixEntry(n, nrows * i, ncols * i)

        self.recalcMatrixSize()
        print(self.rows(), self.cols())


class GeostatisticConstraintsMatrix(pgcore.MatrixBase):
    """Geostatistic constraints matrix

    Uses geostatistical operators described by Jordi et al. (2018),
    however corrects for the remaining non-smooth (damping) part by
    correcting for the spur of the inverse root matrix.

    Jordi, C., Doetsch, J., Günther, T., Schmelzbach, C. & Robertsson, J.O.A.
    (2018): Geostatistical regularisation operators for geophysical inverse
    problems on irregular meshes. Geoph. J. Int. 213, 1374-1386,
    doi:10.1093/gji/ggy055.
    """
    def __init__(self, CM=None, mesh=None, **kwargs):
        """Initialize by computing the covariance matrix & its inverse root.

        Parameters
        ----------
        CM : pg.Matrix or pg.SparseMapMatrix
            covariance matrix, if not given, use mesh and I
        mesh : pg.Mesh
            mesh of which the cell midpoints are used for covariance
        I : float | iterable of floats
            axis correlation length (isotropic) or lengths (anisotropic)
        dip : float [0]
            angle of main axis corresponding to I[0] (2D) or I[0]&I[1] (3D)
        strike : float [0]
            angle of main axis corresponding to I[0] versus I[1] (3D)
        withRef : bool [False]
            neglect spur (reference model effect) that is otherwise corrected
        """
        from pygimli.utils.geostatistics import covarianceMatrix

        super().__init__()
        if isinstance(CM, pgcore.Mesh):
            CM = covarianceMatrix(CM, **kwargs)
        if CM is None:
            CM = covarianceMatrix(mesh, **kwargs)

        self.nModel = CM.shape[0]
        self.CM05 = Cm05Matrix(CM)
        self.spur = self.CM05 * pgcore.RVector(self.nModel, 1.0)
        if kwargs.pop('withRef', False):
            self.spur *= 0.0

    def mult(self, x):
        return self.CM05.mult(x) - self.spur * x

    def transMult(self, x):
        return self.CM05.transMult(x) - self.spur * x

    def cols(self):
        return self.nModel

    def rows(self):
        return self.nModel
