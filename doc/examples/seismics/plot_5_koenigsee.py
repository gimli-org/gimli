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

################################################################################
# We import pyGIMLi and the refraction manager.

import pygimli as pg
from pygimli.physics import TravelTimeManager

################################################################################
# The helper function `pg.getExampleFile` downloads the data set and saves it
# into a temporary location.

data = pg.getExampleFile("traveltime/koenigsee.sgt", load=True, verbose=True)

################################################################################
# We initialize the refraction manager.
mgr = TravelTimeManager()

################################################################################
# Let's have a look at the data in the form of traveltime curves and apparent
# velocity images.
mgr.showData(data)  # show first arrivals as curves (done later with response)
#TODO mgr.showVA(data)  # show data as apparent velocity image

################################################################################
# Finally, we call the `invert` method and plot the result.The mesh is created
# based on the sensor positions on-the-fly. Yes, it is really as simple as that.

mesh = pg.meshtools.createParaMesh(data.sensors(), 
                                   boundary=0, paraMaxCellSize=1.0)
mgr.invert(data, mesh=mesh, secNodes=3, 
           zWeight=0.2, vTop=500, vBottom=5000, 
           verbose=1)
ax, cbar = mgr.showResult()
mgr.showRayPaths(ax=ax, color="0.8", lw=0.3, alpha=0.3)

mgr.showResultAndFit()

################################################################################
# You can play around with the gradient starting model (`vtop` and `vbottom`
# arguments) and the regularization strength `lam`. You can also customize the
# mesh.

pg.wait()