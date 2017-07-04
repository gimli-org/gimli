#!/usr/bin/env python
# encoding: utf-8

r"""
Traveltime in 2D
----------------
"""


import numpy as np

import pygimli as pg

import pygimli.meshtools as mt
from pygimli.physics.traveltime import Refraction

# Create geometry definition for the modelling domain
world = mt.createWorld(start=[-20, 0], end=[20, -16], layers=[-2, -8],
                       worldMarker=False)

# Create a mesh from the geometry definition
mesh = mt.createMesh(world, quality=33, area=0.1, smooth=[1, 10])

pg.show(mesh, mesh.cellMarker(), label='marker')

scheme = pg.physics.traveltime.createRAData(np.linspace(-15., 15., 31))

ra = Refraction()

data = ra.simulate(mesh, np.array(mesh.cellMarker())*1000,
                   scheme, noiseLevel=0.01, noiseAbs=1e-5, verbose=1)
ra.showVA(data)
ra.showData(data)

pg.wait()
