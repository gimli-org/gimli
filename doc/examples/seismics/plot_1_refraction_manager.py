"""
Refraction Manager
------------------

The layered example is taken from the Bachelor thesis of Constanze Reinken
(University of Bonn).
"""

import numpy as np

import pygimli as pg
import pygimli.meshtools as mt
from pygimli.meshtools import (createCircle, createMesh, createPolygon,
                               createRectangle, createWorld, mergePLC)
from pygimli.physics.traveltime import Refraction

################################################################################
# We start by creating a three-layered slope.

layer1 = createPolygon([[0.0, 137], [117.5, 164], [117.5, 162], [0.0, 135]],
                       isClosed=True, marker=1)
layer2 = createPolygon([[0.0, 126], [0.0, 135], [117.5, 162], [117.5, 153]],
                       isClosed=True, marker=2)
layer3 = createPolygon([[0.0, 110], [0.0, 126], [117.5, 153], [117.5, 110]],
                       isClosed=True, marker=3)

slope = (164 - 137) / 117.5

geom = mt.mergePLC([layer1, layer2, layer3])
mesh = mt.createMesh(geom, quality=33, area=0.8, smooth=[1, 10])
pg.show(mesh)

################################################################################
# Next we define geophone positions and a measurement scheme, which consists of
# shot and receiver indices.

numberGeophones = 48
sensors = np.linspace(0., 117.5, numberGeophones)
scheme = pg.physics.traveltime.createRAData(sensors)

# Adapt sensor positions to slope
pos = np.array(scheme.sensorPositions())
for x in pos[:, 0]:
    i = np.where(pos[:, 0] == x)
    new_y = x * slope + 137
    pos[i, 1] = new_y

scheme.setSensorPositions(pos)

################################################################################
# Now we initialize the refraction manager and asssign P-wave velocities to the
# laylers.

ra = Refraction()

vp = np.array(mesh.cellMarkers())
vp[vp == 1] = 250
vp[vp == 2] = 500
vp[vp == 3] = 1300

################################################################################
# We use this model to create noisified synthetic data and look at the
# traveltime curves.

data = ra.simulate(mesh, 1.0 / vp, scheme, noiseLevel=0.01, noiseAbs=0,
                   verbose=True)
ra.showData(data)

################################################################################
# And finally to invert the synthetic data on a new mesh without a priori
# information on the layered structure.

ra = Refraction(data)
ra.createMesh(depth=30., paraMaxCellSize=5.0)
inv = ra.invert()
pg.show(ra.mesh, inv, label="Velocity [m/s]", cMin=250, cMax=1300)
