#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Field data inversion ("Koenigsee")
----------------------------------

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
from pygimli.physics import Refraction

################################################################################
# The helper function `pg.getExampleFile` downloads the data set and saves it
# into a temporary location.

filename = pg.getExampleFile("traveltime/koenigsee.sgt")

################################################################################
# We initialize an instance of the refraction manager with the filename.

ra = Refraction(filename)
print(ra)

################################################################################
# Let's have a look at the data in the form of traveltime curves and apparent
# velocity images.

ra.showData()  # show first arrivals as curves (done later with response)
ra.showVA()  # show data as apparent velocity image

################################################################################
# Finally, we call the `invert` method and plot the result.The mesh is created
# based on the sensor positions on-the-fly. Yes, it is really as simple as that.

ra.invert(zWeight=0.2)
ra.showResult()

################################################################################
# You can play around with the gradient starting model (`vtop` and `vbottom`
# arguments) and the regularization strength `lam`. You can also customize the
# mesh by calling `ra.createMesh()` with options of your choice prior to the
# inversion call.
