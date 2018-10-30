#!/usr/bin/env python
# encoding: utf-8
# sphinx_gallery_thumbnail_number = 2
"""
Refraction Manager
------------------

This example shows how to use the Refraction manager to generate the response
of a three-layered sloping model and to invert the synthetic noisified data."""

import numpy as np

import pygimli as pg
import pygimli.meshtools as mt
from pygimli.physics import Refraction

###############################################################################
# We start by creating a three-layered slope (The model is taken from the BSc
# thesis of Constanze Reinken (University of Bonn).
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

ax, _ = pg.show(mesh, vp, colorBar=True, logScale=False, label='v in m/s')
ax.plot(pos[:, 0], pos[:, 1], 'w+')

###############################################################################
# We use this model to create noisified synthetic data and look at the
# traveltime curves showing three different slopes.

data = ra.simulate(mesh, 1.0 / vp, scheme, noiseLevel=0.001, noiseAbs=0.001,
                   verbose=True)
# ra.showData(data)  # can be used to show the data.

###############################################################################
# And invert the synthetic data on an independent mesh without # information on
# the layered structure we create a new instance of the class using the data.
# Instead of the data, a file name can be given. This is probably where most
# users with data start. See refraction class for supported formats.

ra = Refraction(data)
ra.showData()
ra.createMesh(depth=30., paraMaxCellSize=5.0, secNodes=1)  # can be omitted
# ra.fop.createJacobian(np.ones(ra.fop.regionManager().parameterCount()))
# exit()
vest = ra.invert()  # estimated velocity distribution

###############################################################################
# The method showResult is used to plot the result. Note that only covered
# cells are shown by default. For comparison we plot the geometry on top.

ax, cb = ra.showResult(cMin=min(vp), cMax=max(vp), logScale=False)
pg.show(geom, ax=ax, fillRegion=False, regionMarker=False)  # lines on top

###############################################################################
# Note that internally the following is called
#
# .. code-block:: python
#
#    ax, _ = pg.show(ra.mesh, vest, label="Velocity [m/s]", **kwargs)
#

###############################################################################
# Another useful tool is to show the model along with its respone on the data_

ra.showResultAndFit()

###############################################################################
# Takeaway message:
# A default data inversion with checking of the data consists of only few lines
# (Everthing else can be looked at by introspecting the Refraction manager)
#
# .. code-block:: python
#
#    from pygimli.physics import Refraction
#    ra = Refraction(filename)
#    ra.invert()
#    ra.showResultAndFit()

pg.wait()
