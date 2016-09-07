# -*- coding: utf-8 -*-

"""
Some matrix specialization.
"""

from pygimli.core import _pygimli_ as pg


class LMultRMatrix(pg.MatrixBase):
    """ matrix consisting of actual RMatrix and lef-side vector"""
    def __init__(self, A, left, verbose=False):
        """ constructor saving matrix and vector """
        if A.rows() != len(left):
            raise Exception("Matrix columns do not fit vector length!")

        self.A = A
        self.left = left
        super().__init__(verbose)  # only in Python 3
#        pg.MatrixBase.__init__(self)  # the Python 2 variant

    def rows(self):
        """ return number of rows (using underlying matrix) """
        return self.A.rows()

    def cols(self):
        """ return number of columns (using underlying matrix) """
        return self.A.cols()

    def mult(self, x):
        """ multiplication from right-hand-side (dot product) """
        return self.A.mult(x) * self.left

    def transMult(self, x):
        """ multiplication from right-hand-side (dot product) """
        return self.A.transMult(x * self.left)


class RMultRMatrix(pg.MatrixBase):
    """ Matrix A to be multiplied by a right hand side vector r. """
    def __init__(self, A, r=None):
        super().__init__()
        self.A = A

        if r is None:
            self.r = pg.RVector(A.cols(), 1.0)
        else:
            self.r = r

    def mult(self, x):
        """ return M*x = A*(r*x) """
        return self.A.mult(x*self.r)

    def transMult(self, x):
        """ return (A.T*x)*r """
        return self.A.transMult(x) * self.r

    def cols(self):
        """ number of columns """
        return self.A.cols()

    def rows(self):
        """ number of rows """
        return self.A.rows()


class LRMultRMatrix(pg.MatrixBase):
    """ matrix consisting of actual RMatrix and lef-side vector"""
    def __init__(self, A, left, right, verbose=False):
        """ constructor saving matrix and vector """
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
        """ return number of rows (using underlying matrix) """
        return self.A.rows()

    def cols(self):
        """ return number of columns (using underlying matrix) """
        return self.A.cols()

    def mult(self, x):
        """ multiplication from right-hand-side (dot product) """
        return self.A.mult(x * self.right) * self.left

    def transMult(self, x):
        """ multiplication from right-hand-side (dot product) """
        return self.A.transMult(x * self.left) * self.right
