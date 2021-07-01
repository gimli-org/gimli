#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
ERT field data with topography
------------------------------

Simple example of data measured over a slagdump demonstrating:

- 2D inversion with topography
- geometric factor generation
- topography effect
"""
# sphinx_gallery_thumbnail_number = 3
import pygimli as pg

from pygimli.physics.ert import ERTManager, createGeometricFactors

###############################################################################
# Get some example data with topography
#
data = pg.getExampleFile('ert/slagdump.ohm', load=True, verbose=True)

###############################################################################
# The data file does not contain geometric factors (token field 'k'),
# so we create them based on the given topography.
data['k'] = createGeometricFactors(data, numerical=True)

###############################################################################
# We initialize the ERTManager for further steps and eventually inversion.
ert = ERTManager(sr=False, useBert=True, verbose=True, debug=False)

###############################################################################
# It might be interesting to see the topography effect, i.e the ratio between
# the numerically computed geometry factor and the analytical formula
k0 = createGeometricFactors(data)
ert.showData(data, vals=k0/data['k'], label='Topography effect')

###############################################################################
# The data container has no apparent resistivities (token field 'rhoa') yet.
# We can let the Manager fix this later for us (as we now have the 'k' field),
# or we do it manually.
ert.checkData(data)
print(data)

###############################################################################
# The data container does not necessarily contain data errors data errors
# (token field 'err'), requiring us to enter data errors. We can let the
# manager guess some defaults for us automaticly or set them manually
data['err'] = ert.estimateError(data, absoluteError=0.001, relativeError=0.03)
# or manually:
# data['err'] = data_errors  # somehow

###############################################################################
# Now the data have all necessary fields ('rhoa', 'err' and 'k') so we can run
# the inversion. The inversion mesh will be created with some optional values
# for the parametric mesh generation.
#
mod = ert.invert(data, lam=10,
                 paraDX=0.3, paraMaxCellSize=10, paraDepth=20, quality=33.6)

###############################################################################
# We can view the resulting model in the usual way.
ert.showResultAndFit()
# np.testing.assert_approx_equal(ert.inv.chi2(), 1.10883, significant=3)

###############################################################################
# Or just plot the model only.
ert.showModel(mod)
