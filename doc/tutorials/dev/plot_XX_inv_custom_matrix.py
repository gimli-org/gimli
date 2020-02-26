#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Author: carsten-forty2

"""
import pygimli as pg


class MyMatrix(pg.MatrixBase):
    def __init__(self, verbose=False):
        super().__init__(verbose)

    def rows(self):
        """ Return amount of rows """
        print('musthave rows')
        return 0

    def cols(self):
        """ Return amount of columns """
        print('musthave cols')
        return 0

    def save(self, name):
        """ Save this matrix (optional) but mandatory if used inf inversion
        that is called with dosave=True. """
        print('musthave save if inversion is called with dosave ')
        return 0


class TestModelling(pg.core.ModellingBase):
    def __init__(self, verbose):
        super().__init__(verbose)

        self._J = MyMatrix()
        self.setJacobian(self._J)

    def createJacobian(self, model):
        print('Create Jacobian')

    def response(self, par):
        print('Create response')
        return [0.0, 0.0]

fop = TestModelling(verbose=True)
fop.setStartModel([0, 0])

inv = pg.Inversion(verbose=True, dosave=False)

#fixme!!!!!!!! .. segfault if one of the following is unset
inv.setForwardOperator(fop)
inv.setData([0.0, 0.0])

inv.start()
