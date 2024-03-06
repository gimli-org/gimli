#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
ERT field data with topography
==============================

Simple example of data measured over a slagdump demonstrating:

- 2D inversion with topography
- geometric factor generation
- topography effect
"""
# sphinx_gallery_thumbnail_number = 5
import pygimli as pg
from pygimli.physics import ert

###############################################################################
# Get some example data with topography, typically by a call like
# data = ert.load("filename.dat")
# that supports various file formats
data = pg.getExampleFile('ert/slagdump.ohm', load=True, verbose=True)
print(data)

###############################################################################
# The data file does not contain geometric factors (token field 'k'),
# so we create them based on the given topography.
data['k'] = ert.createGeometricFactors(data, numerical=True)

###############################################################################
# We initialize the ERTManager for further steps and eventually inversion.
mgr = ert.ERTManager(sr=False)

###############################################################################
# It might be interesting to see the topography effect, i.e the ratio between
# the numerically computed geometry factor and the analytical formula
k0 = ert.createGeometricFactors(data)
ert.showData(data, vals=k0/data['k'], label='Topography effect')

###############################################################################
# The data container has no apparent resistivities (token field 'rhoa') yet.
# We can let the Manager fix this later for us (as we now have the 'k' field),
# or we do it manually.
mgr.checkData(data)
print(data)

###############################################################################
# The data container does not necessarily contain data errors data errors
# (token field 'err'), requiring us to enter data errors. We can let the
# manager guess some defaults for us automaticly or set them manually
data['err'] = ert.estimateError(data, relativeError=0.03, absoluteUError=5e-5)
# or manually:
# data['err'] = data_errors  # somehow
ert.show(data, data['err']*100)

###############################################################################
# Now the data have all necessary fields ('rhoa', 'err' and 'k') so we can run
# the inversion. The inversion mesh will be created with some optional values
# for the parametric mesh generation.
#
mod = mgr.invert(data, lam=10, verbose=True,
                 paraDX=0.3, paraMaxCellSize=10, paraDepth=20, quality=33.6)

mgr.showResult()

###############################################################################
# We can view the resulting model in the usual way.
mgr.showResultAndFit()
# np.testing.assert_approx_equal(ert.inv.chi2(), 1.10883, significant=3)

###############################################################################
# Or just plot the model only using your own options.
mgr.showResult(mod, cMin=5, cMax=30, cMap="Spectral_r", logScale=True)
