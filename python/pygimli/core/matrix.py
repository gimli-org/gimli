# -*- coding: utf-8 -*-

"""
Some matrix specialization.
"""

from pygimli.core import _pygimli_ as pg

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