#!/usr/bin/env python
# encoding: utf-8
"""
2D crosshole ERT inversion
--------------------------

Inversion of 2D crosshole field data.
"""

# %%%
# We import the used libraries pygimli, meshtools the ERT module and a function
# for displaying a data container.
#

import numpy as np
import pygimli as pg
import pygimli.meshtools as mt
from pygimli.physics import ert
from pygimli.viewer.mpl import showDataContainerAsMatrix
from pygimli.viewer import pv

# %%%
# We load the data file from the example repository. It represents a crosshole
# data set published by Kuras et al. (2009) in the frame of the ALERT project.
#

data = pg.getExampleData("ert/crosshole3d.dat")
print(data)

# %%%
# There are 36 electrodes, each 9 in four boreholes, and 1256 data that are
# resistances only. Therefore we first compute the geometric factors and then
# the apparent resistivities, of which we plot the last few values.
#

data["k"] = ert.geometricFactors(data, dim=3)
data["rhoa"] = data["r"] * data["k"]

# %%%
# We plot the data in form of a crossplot between the A-B and M-N electrode
# combinations.
#

ab = data["a"] * 100 + data["b"]
mn = data["m"] * 100 + data["n"]
showDataContainerAsMatrix(data, ab, mn, "rhoa", cMap="Spectral_r")

# %%%
# We determine the x and z positions and create a regular grid with a spacing
# of 5xm that contains the electrodes as nodes. This is not necessary but
# improves quality of the forward response. Around the boreholes there is 0.5m
# space and all mesh cells have the marker 2.
#


# %%%
# We first extract the borehole locations, i.e. the x and y positions of the
# electrodes. From these we create a rectangle with 40% boundary and marker 2
# and add the borehole positions to it.
#

elPosXY = np.unique(np.column_stack([pg.x(data), pg.y(data)]), axis=0)
rect = mt.createRectangle(pnts=elPosXY, minBBOffset=1.4, marker=2)
for elpos in elPosXY:
    rect.createNode(*elpos, 0)

ax, cb = pg.show(rect)
ax.plot(*elPosXY.T, "mx")

# %%%
# From this PLC, we create a mesh using a maximum cell size.
# We add an outer (modelling) boundary.
#

bnd = 5
rectMesh = mt.createMesh(rect, quality=34.5, area=.2)
mesh2d = mt.appendTriangleBoundary(
    rectMesh, boundary=bnd, isSubSurface=False, marker=1)
ax, cb = pg.show(mesh2d, markers=True, showMesh=True)
ax.plot(*elPosXY.T, "mx")

# %%%
# We create a vertical discretization vector with dense spacing in the range of
# the electrodes and a coarser discretization above and below. From the 2d mesh
# and the z vector we create a 3D triangular prism mesh that obtains the marker
# of the 2D mesh. Additionally, we set all cells above or below to marker 1
# which is by default the background region.
#

dTop, dBot = 3.5, 10.7
dzIn, dzOut = 0.3, 0.7
zTop = -np.arange(0, dTop, dzOut)  # the upper layer
zMid = -np.arange(dTop, dBot, dzIn)  # the middle
zBot = -np.arange(dBot, dBot+bnd+.1, dzOut)  # the lower layer
zVec = np.concatenate([zTop, zMid, zBot])
print(zVec)
mesh = mt.createMesh3D(mesh2d, zVec, pg.core.MARKER_BOUND_HOMOGEN_NEUMANN,
                       pg.core.MARKER_BOUND_MIXED)

print(mesh)
for c in mesh.cells():
    cd = -c.center().z()  # center depth
    if cd < dTop or cd > dBot:
        c.setMarker(1)

mesh.exportVTK("mesh.vtk")
mesh.exportBoundaryVTU("mesh.vtu")

# %%%
# We estimate an error using default values, i.e. 3% relative error and an
# absolute error of 100uV at an assumed current of 100mA which is almost zero.
# Inversion is run with less weight into the vertical direction.
#

data["err"] = ert.estimateError(data)
mgr = ert.Manager(data)
mgr.invert(mesh=mesh, zWeight=0.4, verbose=True)

# %%%
# We visualize the result
#

pd = mgr.paraDomain
pd["res"] = mgr.model
pd.exportVTK("result.vtk")
pl, _ = pg.show(pd, label="res", style="surface", cMap="Spectral_r",  # hold=1
                filter={"slice": dict(normal=[-1, -1, 0], origin=[2, 2, -6])})
pl.camera_position = "yz"
pl.camera.azimuth = 20
pl.camera.elevation = 20
pl.camera.zoom(1.2)
pl.show()

# %%%
# References
# ----------
#
