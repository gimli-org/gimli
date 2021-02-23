#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
.. _ex:koenigsee:

Field data inversion ("Koenigsee")
==================================

This minimalistic example shows how the Refraction Manager can be used to invert
a field data set. Here, we consider the Koenigsee data set, which represents
classical refraction seismics data set with slightly heterogeneous overburden
and some high-velocity bedrock. The data file can be found in the `pyGIMLi
example data repository
<https://github.com/gimli-org/example-data/blob/master/traveltime/koenigsee.sgt>`_.
"""
# sphinx_gallery_thumbnail_number = 2

################################################################################
# We import pyGIMLi and the refraction manager.

import pygimli as pg
from pygimli.physics import TravelTimeManager

################################################################################
# The helper function `pg.getExampleFile` downloads the data set and saves it
# into a temporary location. Printing the data reveals that there are 714 data
# points using 63 sensors (shots and geophones) with the data columns s (shot),
# g (geophone), and t (traveltime). By default, there is also a validity flag.

data = pg.getExampleFile("traveltime/koenigsee.sgt", load=True, verbose=True)
print(data)

################################################################################
# Let's have a look at the data in the form of traveltime curves.

fig, ax = pg.plt.subplots()
pg.physics.traveltime.drawFirstPicks(ax, data)

################################################################################
# We initialize the refraction manager.
mgr = TravelTimeManager()

# Alternatively, one can plot a matrix plot of apparent velocities which is the
# more general function also making sense for crosshole data.
ax, cbar = mgr.showData(data)

################################################################################
# Finally, we call the `invert` method and plot the result.The mesh is created
# based on the sensor positions on-the-fly.

mgr.invert(data, secNodes=3, paraMaxCellSize=5.0,
           zWeight=0.2, vTop=500, vBottom=5000,
           verbose=1)

ax, cbar = mgr.showResult(logScale=True)
mgr.drawRayPaths(ax=ax, color="w", lw=0.3, alpha=0.5)

################################################################################
# Show result and fit of measured data and model response. You may want to save your results too.
fig = mgr.showResultAndFit()
mgr.saveResult()
################################################################################
# You can play around with the gradient starting model (`vTop` and `vBottom`
# arguments) and the regularization strength `lam`. You can also customize the
# mesh.
