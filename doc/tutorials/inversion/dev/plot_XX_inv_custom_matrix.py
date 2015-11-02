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
        print('musthave rows')
        return 0

    def cols(self):
        print('musthave cols')
        return 0

    def save(self, name):
        print('musthave save')
        return 0
    
        
class TestModelling(pg.ModellingBase):
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

inv = pg.RInversion(verbose=True, dosave=False)
inv.setForwardOperator(fop)
inv.setData([0.0,0.0])

inv.start()

