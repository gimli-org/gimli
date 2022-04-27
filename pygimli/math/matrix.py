# -*- coding: utf-8 -*-
"""Some special and usefull matrices."""

import numpy as np

import pygimli as pg
import pygimli.core as pgcore


from pygimli.core import (CMatrix, CSparseMapMatrix, CSparseMatrix,
                          RSparseMapMatrix, RSparseMatrix, ElementMatrix,
                          IVector, MatrixBase, R3Vector, RVector)

IdentityMatrix = pgcore.IdentityMatrix
Matrix = pgcore.RMatrix
SparseMatrix = pgcore.RSparseMatrix
SparseMapMatrix = pgcore.RSparseMapMatrix
BlockMatrix = pgcore.RBlockMatrix



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
        """Initialize matrix by matrix and vectors.

        Parameters
        ----------
        A : pg.core.MatrixBase derived
            main matrix
        left, right : pg.Vector
            left and right side vectors to weight individual rows & columns
        """
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
    """Matrix by addition of two matrices implicitly.

        The matrix holds two matrices and distributes multiplication in 2 parts

        M = A + B
        M * x = A * x + B * x
        M.T * y = B^T * y + A.T^*y
    """

    def __init__(self, A, B):
        """Initialize matrix.

        Parameters
        ----------
        A, B : any Matrix derived from pg.core.MatrixBase
            submatrices to be multiplied (A/B.cols() and rows() need to match)
        """
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
    """Matrix by multipication of two matrices implicitly.

        The matrix holds two matrices and distributes multiplication in 2 parts

        M = A * B
        M * x = A * (B*x)
        M.T * y = B^T * (A.T^*y) = (y^T * A * B)^T
    """

    def __init__(self, A, B):
        """Initialize matrix.

        Parameters
        ----------
        A, B : any Matrix derived from pg.core.MatrixBase
            submatrices to be multiplied (A.cols() must equal B.rows())
        """
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
        """Initialize matrix

        Parameters
        ----------
        d : Vector
            vector holding diagonal elements
        """
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

@pg.cache
def createCm05(A):
    """Globally cached helper function to create Cm05Matrix."""
    return pg.matrix.Cm05Matrix(A, verbose=True)

class Cm05Matrix(pgcore.MatrixBase):
    """Matrix implicitly representing the inverse square-root."""

    def __init__(self, A, verbose=False):
        """Constructor saving matrix and vector.

        Parameters
        ----------
        A : ndarray
            numpy type (full) matrix
        """
        super().__init__(verbose)  # only in Python 3
        self._mul = None

        if isinstance(A, str):
            self.load(A)
        else:
            from scipy.linalg import eigh  # , get_blas_funcs

            if A.shape[0] != A.shape[1]:  # rows/cols for pgcore matrix
                raise Exception("Matrix must by square (and symmetric)!")

            if verbose:
                t = pg.tic(key='init cm05')
            self.ew, self.EV = eigh(A)

            if verbose:
                pg.info('(C) Time for eigenvalue decomposition: {:.1f}s'.format(
                    pg.dur(key='init cm05')))

            #self.A = A

    @property
    def mul(self):
        if self._mul is None:
            self._mul = np.sqrt(1./self.ew)
        return self._mul

    def save(self, fileName):
        """Save the content of this matrix. Used for caching until pickling is possible for this class"""
        np.save(fileName, dict(ew=self.ew, EV=self.EV), allow_pickle=True)

    def load(self, fileName):
        """Load the content of this matrix. Used for caching until pickling is possible for this class"""
        d = np.load(fileName + '.npy', allow_pickle=True).tolist()
        self.ew = d['ew']
        self.EV = d['EV']

    def rows(self):
        """Return number of rows (using underlying matrix)."""
        return len(self.ew)

    def cols(self):
        """Return number of columns (using underlying matrix)."""
        return self.row()

    def mult(self, x):
        """Multiplication from right-hand side (dot product)."""
        part1 = (np.dot(np.transpose(x), self.EV).T*self.mul).reshape(-1, 1)
        return self.EV.dot(part1).reshape(-1,)
#        return self.EV.dot((x.T.dot(self.EV)*self.mul).T)

    def transMult(self, x):
        """Multiplication from right-hand side (dot product)."""
        return self.mult(x)  # matrix is symmetric by definition


class RepeatVMatrix(pgcore.BlockMatrix):
    """Matrix holding a base matrix N times vertically."""
    def __init__(self, A, num):
        """Initialize."""
        super().__init__()
        self.A_ = A
        self.Aid = self.addMatrix(self.A_)
        nr = 0
        for i in range(num):
            self.addMatrixEntry(self.Aid, nr, 0)
            nr += A.rows()

        self.recalcMatrixSize()


class RepeatHMatrix(pgcore.BlockMatrix):
    """Matrix holding a base matrix N times vertically."""
    def __init__(self, A, num):
        """Initialize."""
        super().__init__()
        self.A_ = A
        self.Aid = self.addMatrix(self.A_)
        nc = 0
        for i in range(num):
            self.addMatrixEntry(self.Aid, 0, nc)
            nc += A.cols()

        self.recalcMatrixSize()


class RepeatDMatrix(pgcore.BlockMatrix):
    """Matrix holding a base matrix N times vertically."""
    def __init__(self, A, num):
        """Initialize."""
        super().__init__()
        self.A_ = A
        self.Aid = self.addMatrix(self.A_)
        nc = 0
        nr = 0
        for i in range(num):
            self.addMatrixEntry(self.Aid, nr, nc)
            nc += A.cols()
            nr += A.rows()

        self.recalcMatrixSize()


class NDMatrix(pgcore.BlockMatrix):
    """Diagonal block (block-Jacobi) matrix derived from pg.matrix.BlockMatrix.

    (to be moved to a better place at a later stage)
    """

    def __init__(self, num, nrows, ncols):
        super().__init__()  # call inherited init function
        self.Ji = []  # list of individual block matrices
        for i in range(num):
            self.Ji.append(pgcore.Matrix())
            self.Ji[-1].resize(nrows, ncols)
            n = self.addMatrix(self.Ji[-1])
            self.addMatrixEntry(n, nrows * i, ncols * i)

        self.recalcMatrixSize()


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
        super().__init__(kwargs.pop('verbose', False))
        self.withRef = kwargs.pop('withRef', False)
        self._spur = None

        if isinstance(CM, str):
            self.load(CM)
        else:
            from pygimli.utils.geostatistics import covarianceMatrix

            if isinstance(CM, pgcore.Mesh):
                CM = covarianceMatrix(CM, **kwargs)

            if CM is None and mesh is not None:
                CM = covarianceMatrix(mesh, **kwargs)
            else:
                pg.critical('Give either CM or mesh')

            self.Cm05 = createCm05(CM)

    @property
    def spur(self):
        if self._spur is None:
            if self.withRef is True:
                self._spur = np.zeros(self.rows())
            else:
                self._spur = self.Cm05 * pgcore.RVector(self.rows(), 1.0)
        return self._spur

    @property
    def nModel(self):
        try:
            return self.Cm05.size()
        except Exception as e:
            return 0

    def save(self, fileName):
        """Save the content of this matrix. Used for caching until pickling is possible for this class
        """
        self.Cm05.save(fileName + '-Cm05')
        np.save(fileName, dict(verbose=self.verbose(),
                               withRef=self.withRef,
                               Cm05=fileName +'-Cm05'),
                        allow_pickle=True)


    def load(self, fileName):
        """Load the content of this matrix. Used for caching until pickling is possible for this class
        """
        d = np.load(fileName + '.npy', allow_pickle=True).tolist()
        self.setVerbose(d['verbose'], )
        self.withRef = d['withRef']
        self.Cm05 = Cm05Matrix(d['Cm05'])


    def mult(self, x):
        return self.Cm05.mult(x) - self.spur * x

    def transMult(self, x):
        return self.Cm05.transMult(x) - self.spur * x

    def cols(self):
        return self.nModel

    def rows(self):
        return self.nModel

    def clear(self):
        self.Cm05 = None
        self._spur = None
