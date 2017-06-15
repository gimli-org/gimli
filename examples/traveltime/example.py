#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Example of reading, displaying and inverting refraction data using"""
from pygimli.physics import Refraction

ra = Refraction('koenigsee.sgt')
print(ra)
ra.showData()  # show first arrivals as curves (done later with response)
ra.showVA()  # show data as apparent velocity image
if False:  # user defined mesh
    ra.createMesh()  # depth=<user-defined depth>, quality=<trangle quality>
    ra.showMesh()  # show just the mesh for checking

# ra.estimateError(absoluteError=0.0006, relativeError=0.001)
ra.invert(zWeight=0.2)  #
# options: vtop/vbottom..gradient starting model velocity, zweight=0.2, lam=30
# ra.showResult()  # velocity image, cMin/cMax - color range, logScale=bool
ra.showResultAndFit()
