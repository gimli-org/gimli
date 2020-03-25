#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
Simple example -- slagdump testcase
- 2d with topo
- geometric factor generation
- topography effect
"""
import numpy as np
import pygimli as pg

import pygimli.meshtools as mt
from pygimli.physics.ert import ERTManager, createGeometricFactors

##############################################################################
# Get some example data with topogrpahy
data = pg.getExampleFile('ert/slagdump.ohm', load=True, verbose=True)

##############################################################################
# Initialize the ERTManager
ert = ERTManager(sr=False, useBert=True, verbose=True, debug=False)
                      
##############################################################################
# We need a mesh for every inversion
mesh = mt.createParaMesh(data.sensors(), 
                         paraDX=0.3, paraMaxCellSize=10, paraDepth=20,
                         quality=33.6)

##############################################################################
# The data file does not contain geometric factors so create them for the given 
# topography and show the topography effect
kTopo = createGeometricFactors(data, mesh)
k0 = createGeometricFactors(data)

ert.showData(data, vals=k0/kTopo, label='Topography effect')

data['k'] = kTopo
##############################################################################
# data have no rhoa .. let the manager fix this (will be done automatic)
# ert.dataCheck(data)

##############################################################################
# data have no err .. let the manager guess some (will be done automatic)
# data['err'] = ert.estimateError(data, absoluteError=0.001, relativeError=0.03)

mod = ert.invert(data, mesh=mesh, maxIter=20, lam=10)
ert.showResultAndFit()
np.testing.assert_approx_equal(ert.inv.chi2(), 1.10883, significant=3)

pg.wait()
