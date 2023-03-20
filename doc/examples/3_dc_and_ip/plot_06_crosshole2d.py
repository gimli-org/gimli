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

# %%%
# We load the data file from the example repository. It represents a crosshole
# data set published by Kuras et al. (2009) in the frame of the ALERT project.
#

# data = ert.load("crosshole2d.dat")
data = pg.getExampleData("ert/crosshole2d.dat")
print(data)

# %%%
# There are 144 electrodes, each 16 in nine boreholes, and 1256 data that are
# resistances only. Therefore we first compute the geometric factors and then
# the apparent resistivities, of which we plot the last few values.
#

data["k"] = ert.geometricFactors(data, dim=3)
data["rhoa"] = data["r"] * data["k"]
print(np.sort(data["rhoa"])[-5:])

# %%%
# We see that there is a single value above 500 Ohmm and a few other high ones.
# Therefore we delete the highest value. We plot the data in form of a
# crossplot between the A and M electrodes (with forward and backward sign).
#

data.remove(data["rhoa"] > 200)
m = np.sign(data["m"] - data["n"]) * data["m"]
showDataContainerAsMatrix(data, m, data["a"], "rhoa", cMap="Spectral_r")

# %%%
# We determine the x and z positions and create a regular grid with a spacing
# of 5xm that contains the electrodes as nodes. This is not necessary but
# improves quality of the forward response. Around the boreholes there is 0.5m
# space and all mesh cells have the marker 2.
#

ex = np.unique(pg.x(data))
ez = np.unique(pg.z(data))
dx = 0.05
nb = 8
xmin, xmax = min(ex) - nb*dx, max(ex) + nb*dx
zmin, zmax = min(ez) - nb*dx, 0
x = np.arange(xmin, xmax+.001, dx)
z = np.arange(zmin, zmax+.001, dx)
grid = mt.createGrid(x, z, marker=2)
ax, cb = pg.show(grid)
ax.plot(pg.x(data), pg.z(data), "mx")
print(grid)

# %%%
# In order to ensure correct boundary conditions, we append a triangular mesh
# outside that is automatically treated as background because they have a
# region marker of 1 and so-called world boundary conditions, i.e. Neumann BC
# at the Earth's surface and mixed boundary conditions at the other boundaries.
#

mesh = mt.appendTriangleBoundary(grid, marker=1, boundary=5, worldMarkers=1)
pg.show(mesh, markers=True)

# %%%
# We estimate an error using default values, i.e. 3% relative error and an
# absolute error of 100uV at an assumed current of 100mA which is almost zero.
# Inversion is run with half as much weight into the vertical direction.
#

data["err"] = ert.estimateError(data)
mgr = ert.Manager(data)
mgr.inv.setRegularization(correlationLengths=[1, 0.5])
mgr.invert(mesh=mesh, verbose=True)
mgr.showResult(cMin=15, cMax=200)

# %%%
# References
# ----------
# Kuras, O., Pritchard, J., Meldrum, P. I., Chambers, J. E., Wilkinson, P. B.,
# Ogilvy, R. D.,and Wealthall, G. P. (2009). Monitoring hydraulic processes
# with Automated time-Lapse Electrical Resistivity Tomography (ALERT).
# Compte Rendus Geosciences - Special issue on Hydrogeophysics,
# 341(10-11):868–885.
