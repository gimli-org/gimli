# -*- coding: utf-8 -*-
"""Some matrix specialization."""

import time
from pygimli.core import _pygimli_ as pg
import numpy as np

# make core matrices (now in pg, later pg.core) known here for tab-completion
# BlockMatrix = pg.BlockMatrix
# IdentityMatrix = pg.IdentityMatrix


class MultMatrix(pg.MatrixBase):
    def __init__(self, A, verbose=False):
        self._A = A
        super(MultMatrix, self).__init__(verbose)  # only in Python 3

    def rows(self):
        """Return number of rows (using underlying matrix)."""
        return self._A.rows()

    def cols(self):
        """Return number of columns (using underlying matrix)."""
        return self._A.cols()

    def save(self, filename):
        """So it can be used in inversion with dosave flag"""
        pass


class MultLeftMatrix(MultMatrix):
    """Matrix consisting of actual RMatrix and lef-side vector."""

    def __init__(self, A, left, verbose=False):
        """Constructor saving matrix and vector."""
        if A.rows() != len(left):
            raise Exception("Matrix columns do not fit vector length!")
        super(MultLeftMatrix, self).__init__(A, verbose)  # only in Python 3

        self._l = left

    def mult(self, x):
        """Multiplication from right-hand-side (dot product A*x)."""
        return self._A.mult(x) * self._l

    def transMult(self, x):
        """Multiplication from right-hand-side (dot product A.T * x)"""
        return self._A.transMult(x * self._l)


LMultRMatrix = MultLeftMatrix  # alias for backward compatibility


class MultRightMatrix(MultMatrix):
    """Some Matrix, multiplied with a right hand side vector r."""

    def __init__(self, A, r=None, verbose=False):
        super(MultRightMatrix, self).__init__(A, verbose)

        if r is None:
            self._r = pg.RVector(self.cols(), 1.0)
        else:
            self._r = r

    def mult(self, x):
        """Return M*x = A*(r*x)"""
        return self._A.mult(x * self._r)

    def transMult(self, x):
        """Return M.T*x=(A.T*x)*r"""
        return self._A.transMult(x) * self._r


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

    def mult(self, x):
        """Multiplication from right-hand-side (dot product A*x)."""
        return self._A.mult(x * self._r) * self._l

    def transMult(self, x):
        """Multiplication from right-hand-side (dot product A.T*x)."""
        return self._A.transMult(x * self._l) * self._r


LRMultRMatrix = MultLeftRightMatrix  # alias for backward compatibility


class Add2Matrix(pg.MatrixBase):
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


class Mult2Matrix(pg.MatrixBase):
    """Matrix  by multiplying two matrices."""

    def __init__(self, A, B):
        super().__init__()
        self.A = A
        self.B = B
        assert A.cols() == B.rows()

    def mult(self, x):
        """Return M*x = A*(r*x)"""
        return self.A.mult(self.B.mult(x))

    def transMult(self, x):
        """Return M.T*x=(A.T*x)*r"""
        return self.B.transMult(self.A.transMult(x))

    def cols(self):
        """Number of columns."""
        return self.B.cols()

    def rows(self):
        """Number of rows."""
        return self.A.rows()


class DiagonalMatrix(pg.MatrixBase):
    """Square matrix with a vector on the main diagonal."""

    def __init__(self, d):
        super().__init__()
        self.d = d

    def mult(self, x):
        """Return M*x = r*x (element-wise)"""
        return x * self.d

    def transMult(self, x):
        """Return M.T*x=(A.T*x)*r"""
        return x * self.d

    def cols(self):
        """Number of columns (length of diagonal)."""
        return len(self.d)

    def rows(self):
        """Number of rows (length of diagonal)."""
        return len(self.d)


class Cm05Matrix(pg.MatrixBase):
    """Matrix implicitly representing the inverse square-root."""

    def __init__(self, A, verbose=False):
        """Constructor saving matrix and vector.

        Parameters
        ----------
        A : ndarray
            numpy type (full) matrix
        """
        from scipy.linalg import eigh  # , get_blas_funcs

        if A.shape[0] != A.shape[1]:  # rows/cols for pg matrix
            raise Exception("Matrix must by square (and symmetric)!")

        self.size = A.shape[0]
        t = time.time()
        self.ew, self.EV = eigh(A)
        self.mul = np.sqrt(1./self.ew)
        if verbose:
            pg.info('(C) Time for eigenvalue decomposition:{:.1f} s'.format(
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
