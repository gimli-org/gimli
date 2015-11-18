#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Example of reading, displaying and inverting refraction data """

import pygimli as pg
from pygimli.physics.traveltime import Refraction

ra = Refraction('koenigsee.sgt')
print(ra)
#pg.showLater(True)
ra.showData()  # show first arrivals as curves (done later with response)
ra.createMesh()  # depth=<for user-defined depth>, quality=<trangle quality>
ra.showMesh()  # show just the mesh for checking
ra.estimateError(0.0006)  # absolute error (default=1ms), relativeError
ra.run()  # create inversion class and run actual inversion
# options: vtop/vbottom..gradient starting model velocity, zweight=0.2, lam=30
ra.showData(response=ra.inv.response())  # data as x, response as lines
ra.showResult()  # show velocity image, cMin/cMax - color range, logScale=bool
