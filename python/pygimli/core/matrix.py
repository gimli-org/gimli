# -*- coding: utf-8 -*-
"""Some matrix specialization."""

import time
from pygimli.core import _pygimli_ as pg
import numpy as np


class MultLeftMatrix(pg.MatrixBase):
    """Matrix consisting of actual RMatrix and lef-side vector."""

    def __init__(self, A, left, verbose=False):
        """Constructor saving matrix and vector."""
        if A.rows() != len(left):
            raise Exception("Matrix columns do not fit vector length!")

        self.A = A
        self.left = left
        super().__init__(verbose)  # only in Python 3
#        pg.MatrixBase.__init__(self)  # the Python 2 variant

    def rows(self):
        """Return number of rows (using underlying matrix)."""
        return self.A.rows()

    def cols(self):
        """Return number of columns (using underlying matrix)."""
        return self.A.cols()

    def mult(self, x):
        """Multiplication from right-hand-side (dot product A*x)."""
        return self.A.mult(x) * self.left

    def transMult(self, x):
        """Multiplication from right-hand-side (dot product A.T * x)"""
        return self.A.transMult(x * self.left)

LMultRMatrix = MultLeftMatrix  # alias for backward compatibility


class MultRightMatrix(pg.MatrixBase):
    """Some Matrix, multiplied with a right hand side vector r."""

    def __init__(self, A, r=None):
        super().__init__()
        self.A = A

        if r is None:
            self.r = pg.RVector(A.cols(), 1.0)
        else:
            self.r = r

    def mult(self, x):
        """Return M*x = A*(r*x)"""
        return self.A.mult(x * self.r)

    def transMult(self, x):
        """Return M.T*x=(A.T*x)*r"""
        return self.A.transMult(x) * self.r

    def cols(self):
        """Number of columns."""
        return self.A.cols()

    def rows(self):
        """Number of rows."""
        return self.A.rows()

RMultRMatrix = MultRightMatrix  # alias for backward compatibility


class MultLeftRightMatrix(pg.MatrixBase):
    """Matrix consisting of actual RMatrix and left-hand-side vector."""

    def __init__(self, A, left, right, verbose=False):
        """Constructor saving matrix and vector."""
        if A.cols() != len(right):
            raise Exception("Matrix columns do not fit right vector length!")
        if A.rows() != len(left):
            raise Exception("Matrix rows do not fit left vector length!")

        self.A = A
        self.right = right
        self.left = left
        super().__init__(verbose)  # only in Python 3
#        pg.MatrixBase.__init__(self)  # the Python 2 variant

    def rows(self):
        """Number of rows (using the underlying matrix)."""
        return self.A.rows()

    def cols(self):
        """Number of columns (using the underlying matrix)."""
        return self.A.cols()

    def mult(self, x):
        """Multiplication from right-hand-side (dot product A*x)."""
        return self.A.mult(x * self.right) * self.left

    def transMult(self, x):
        """Multiplication from right-hand-side (dot product A.T*x)."""
        return self.A.transMult(x * self.left) * self.right

LRMultRMatrix = MultLeftRightMatrix  # alias for backward compatibility


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
        elapsed = time.time() - t
        print('(C) Calculation time for eigenvalue decomposition:\n%s sec'
              % elapsed)

        self.A = A
        super().__init__(verbose)  # only in Python 3
#        pg.MatrixBase.__init__(self)  # the Python 2 variant

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
