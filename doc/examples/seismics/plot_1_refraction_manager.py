"""
Refraction Manager
------------------

The layered example is taken from the Bachelor thesis of Constanze Reinken
(University of Bonn).
"""

import numpy as np

import pygimli as pg
import pygimli.meshtools as mt
from pygimli.physics.traveltime import Refraction

###############################################################################
# We start by creating a three-layered slope.

layer1 = mt.createPolygon([[0.0, 137], [117.5, 164], [117.5, 162], [0.0, 135]],
                          isClosed=True, marker=1, area=1)
layer2 = mt.createPolygon([[0.0, 126], [0.0, 135], [117.5, 162], [117.5, 153]],
                          isClosed=True, marker=2)
layer3 = mt.createPolygon([[0.0, 110], [0.0, 126], [117.5, 153], [117.5, 110]],
                          isClosed=True, marker=3)

slope = (164 - 137) / 117.5

geom = mt.mergePLC([layer1, layer2, layer3])
mesh = mt.createMesh(geom, quality=34.3, area=3, smooth=[1, 10])
pg.show(mesh)

###############################################################################
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

###############################################################################
# Now we initialize the refraction manager and asssign P-wave velocities to the
# layers. To this end, we create a map from cell markers 0 through 3 to
# velocities (in m/s) and generate a velocity vector.  To check whether the
# it is correct, we plot it # along with the sensor positions

ra = Refraction()
vp = np.array(mesh.cellMarkers())
vp[vp == 1] = 250
vp[vp == 2] = 500
vp[vp == 3] = 1300
# not really nice, this should work as well
# vpmap = np.array([0, 250, 500, 1300], dtype=np.float)
# vp = vpmap[mesh.cellMarkers()]

ax, _ = pg.show(mesh, vp)
ax.plot(pos[:, 0], pos[:, 1], 'w+')

###############################################################################
# We use this model to create noisified synthetic data and look at the
# traveltime curves showing three different slopes.

data = ra.simulate(mesh, 1.0 / vp, scheme, noiseLevel=0.001, noiseAbs=0.001,
                   verbose=True)
ra.showData(data)

###############################################################################
# And finally to invert the synthetic data on an independent mesh without prior
# information on the layered structure.

ra = Refraction(data)
ra.createMesh(depth=30., paraMaxCellSize=5.0)
vest = ra.invert()

###############################################################################
# The method showResult is used to plot the result. Note that only covered
# cells are shown by default. For comparison we plot the geometry on top.
ax, cb = ra.showResult(cMin=min(vp), cMax=max(vp), logScale=False)
pg.show(geom, ax=ax, fillRegion=False, regionMarker=False)
# Note that we could also plot the mesh by hand using
# ax, _ = pg.show(ra.mesh, vest, label="Velocity [m/s]",
#                 cMin=min(vp), cMax=max(vp), logScale=False)
