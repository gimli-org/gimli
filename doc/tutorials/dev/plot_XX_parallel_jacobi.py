#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Author: carsten-forty2

"""
import pygimli as pg
import numpy as np

import time

class TestModelling(pg.core.ModellingBase):
    def __init__(self, nPars, verbose):
        pg.core.ModellingBase.__init__(self, verbose)
        self.regionManager().setParameterCount(nPars)

    def response(self, par):
        print('Create response', str(par))
        time.sleep(1)
        return par * 3.0
    
    def response_mt(self, par, i=0):
        """
            this need to be implemented as read only function!
            don't use self.mapModel until this is changed into a read only 
            version
        """
        print('Create response_mt(' +str(i) +')', str(par))
        time.sleep(1)
        return par * 2.0

nPars = 4
m=pg.Vector(nPars, 1)

fop = TestModelling(nPars, verbose=True)

fop.setMultiThreadJacobian(1)
pg.tic()
fop.createJacobian(m)
pg.toc()
print(fop.jacobian())

fop.setMultiThreadJacobian(4)
pg.tic()
fop.createJacobian(m)
pg.toc()
print(fop.jacobian())