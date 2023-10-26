#!/usr/bin/env python
# encoding: utf-8
"""
# Timelapse ERT
Since pyGIMLi we have a class dedicated to timelapse ERT processing & inversion.
The inversion supports different schemes from simple individual over constrained
inversion to fully coupled ("4D") inversion using pyGIMLi's `MultiFrameModelling`
modelling framework.
Additionally, we created a github repository
https://github.com/gimli-org/timelapseERT that holds published data and scripts
demonstrating how to achieve the published results, according to the FAIR data
standards.

This notebook is a simplistic model for a synthetic case.
"""
# sphinx_gallery_thumbnail_number = 5

# %%%
# We import the used libraries pygimli, meshtools the ERT module.
#

import numpy as np
import matplotlib.pyplot as plt
import pygimli as pg
import pygimli.meshtools as mt
from pygimli.physics import ert

# %%%
# We create a data container with a dipole-dipole array.
#

scheme = ert.createData(elecs=42, spacing=0.5, schemeName='dd',
                        maxSeparation=14)
print(scheme)

# %%%
# Our subsurface is a three-layer model with an aquifer in the middle, into
# which synthetic tracer is injected that is moving to the right.
#

world = mt.createWorld(start=[-50, 0], end=[100, -50], boundary=1,
                       layers=[-1, -5], worldMarker=True)
for pos in scheme.sensorPositions():
    world.createNode(pos, marker=-99)
    world.createNode(pos+pg.RVector3(0, -0.2))

# Create some heterogeneous block
plcs = [world]
pos = [5, 6, 9, 15]
nT = len(pos) - 1  # number of time steps
for i in range(nT):
    block = mt.createRectangle(start=[pos[i], -3], end=[pos[i+1], -2],
                               area=0.1, marker=4+i)
    plcs.append(block)

geom = mt.mergePLC(plcs)
mesh = mt.createMesh(geom, quality=34.4)
print(mesh)
ax, _ = pg.show(mesh, markers=True, boundaryMarkers=False, showMesh=True)
ax.set_xlim(0, 20)
ax.set_ylim(-8, 0)

# %%%
# We associate 100 Ohmm to the first layer, 50 to the second and 20 to the last
# one. At the beginning, the anomalies have the same resistivity as the aquifer.
#

rhomap = [[1, 100.0], [2, 50.0], [3, 20.0], [4, 50.0], [5, 50.0], [6, 50.0]]
noise = dict(noiseLevel=0.01, noiseAbs=0, verbose=False)
mgr = ert.Manager()
data = ert.simulate(mesh=mesh, res=rhomap, scheme=scheme, **noise)
rhoTracer = 10
cDict = dict(colorBar=False, cMin=10, cMax=100, logScale=1, cMap='Spectral_r')
fig, ax = pg.plt.subplots(figsize=(10, 6), ncols=nT+1, nrows=3)
DATA = []
for i in range(nT+1):
    pg.show(mesh, rhomap, ax=ax[0, i], **cDict)
    ax[0, i].set_xlim(0, 20)
    ax[0, i].set_ylim(-8, 0)
    data = ert.simulate(mesh, res=rhomap, scheme=scheme, **noise)
    data.save('data{:d}.dat'.format(i))
    DATA.append(data)
    ert.show(data, ax=ax[1, i], **cDict)
    ratio = data('rhoa') / DATA[0]('rhoa')
    if i > 0:
        ert.show(data, ratio, ax=ax[2, i],
                cMap='bwr', cMin=1/1.5, cMax=1.5, colorBar=False)
    if i < nT:
        rhomap[3+i][1] = rhoTracer

cDict.pop('colorBar')
cDict['label'] = r'$\rho$ [$\Omega$m]'
pg.viewer.mpl.colorbar.createColorBarOnly(ax=ax[2, 0], **cDict)
# ax[2, 0].set_aspect(0.3)

# %%%
# We initialize the `TimelapseERT` class by passing the list of data.
# Otherways are passing
# * a single data file that either holds all timestep or is accompagnied by
# another file with the apparent resistivites (and optionally errors) as matrix
# * a file name with a * in it that points to a number of data (e.g. `bla*.dat`)
# to be read sequentially
# Note that in the latter case or when passing a list of `DataContainerERT` the
# data will be homogenized, i.e. brought to a single `DataContainerERT` and a
# apparent resistivity matrix where non-existing values are masked out.
# In the initialization, one can pass a list of `datatime` objects to specify
# measuring times. Otherwise they are retrieved from the filenames or just set
# to equidistant intervals.
#

tl = ert.TimelapseERT(DATA)
print(tl)

# %%%
# Masking of data can be achieved by `tl.mask()` specifying minimum and maximum
# apparent resistivity or maximum error.
# Additionally, you can filter the data by `tl.filter()` to set maximum geometric
# factor or remove/select timesteps.
# This command can generate multi-page pdfs of the data.
# tl.generateDataPDF(cMin=30, cMax=100)
#

# %%%
# One can do a single timestep inversion using `tl.invert(t=i)`. If this
# argument is omitted, are timesteps are inverted sequentially, always using
# the ERT Manager.
# By default, the first model is used as reference, i.e. the model difference
# is constrained. This can be deactivatey by `isReference=False`.
# The reference model can be moved along with the inversion by by `creep=True`
# so that the difference to the preceding step is constrained.
# One can specify regularization options using a `reg` dictionary and, if
# wanted, a different regularization for the timesteps by `regTL`.
#

tl.invert(zWeight=0.3)
print(tl.chi2s)

# %%%
# After inversion, on can show the models by `showAllModels`. For convenience,
# one can also generate a multi-page pdf by `generateModelPDF`.
#

tl.showAllModels();

# %%%
# Often, one is interested in the changes or differences, which are in the
# usual logarithmic scale the ratios.
#

ax = tl.showAllModels(ratio=True, rMax=3)

# %%%
# As powerful alternative to a sequential inversion, one can invert all
# timesteps together with constraints along the spatial and temporal
# dimensions. For this there is a special call `fullInversion` that might
# take more memory, but is usually not slower than a sequential inversion
# and moreover more robust.
#

tl.fullInversion(zWeight=0.3, lam=100)
ax = tl.showAllModels(ratio=True, rMax=3)

# %%%
# References
# ----------
# To be filled.
