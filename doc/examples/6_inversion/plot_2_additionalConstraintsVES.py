#!/usr/bin/env python
# coding: utf-8
r"""
Incorporating additional constraints as equations
=================================================
Sometimes we have additional information on the model that can be expressed in
form of an equation. Let us consider a simple vertical electrical sounding
(VES) that is inverted for a four-layer case with unknown layer thicknesses.

Assume we know the depth of the layer boundary, e.g. from a borehole, in this
case between the third and fourth layer be :math:`z_2`. So we can formulate this
by an equation :math:`d_1+d_2+d_3=z_2`. This equation can be added to the inverse
problem in the style of Lagrangian, i.e. by an additional regularization term.

We use the LSQR inversion framework as used by Wagner et al. (2019) to
constrain the sum of water, ice, air and rock fraction to be 1.
"""

# sphinx_gallery_thumbnail_number = 2
import numpy as np
import pygimli as pg
from pygimli.frameworks.lsqrinversion import LSQRInversion
from pygimli.physics.ves import VESModelling
from pygimli.viewer.mpl import drawModel1D


# %%%
# We set up a synthetic model of four layers and show the sounding curve
#

nlay = 4  # number of layers
lam = 200.  # (initial) regularization parameter
errPerc = 3.  # relative error of 3 percent
ab2 = np.logspace(-1, 2, 50)  # AB/2 distance (current electrodes)
mn2 = ab2 / 3.  # MN/2 distance (potential electrodes)
f = VESModelling(ab2=ab2, mn2=mn2, nLayers=nlay)
synres = [100., 500., 20., 800.]  # synthetic resistivity
synthk = [0.5, 3.5, 6.]  # synthetic thickness (nlay-th layer is infinite)
rhoa = f(synthk+synres)
rhoa = rhoa * (pg.randn(len(rhoa)) * errPerc / 100. + 1.)
pg.plt.loglog(rhoa, ab2, "x-")
pg.plt.grid(True)

# %%%
# Next, we set up an inversion instance with log transformation on data and
# model side and run the inversion.
#

inv = LSQRInversion(fop=f, verbose=True)
tLog = pg.trans.TransLog()  #
inv.dataTrans = tLog
inv.modelTrans = tLog
startModel = pg.cat(pg.Vector(nlay-1, 8), pg.Vector(nlay, pg.median(rhoa)))
inv.inv.setMarquardtScheme()
error = pg.Vector(len(rhoa), errPerc/100)
model1 = inv.run(rhoa, error, lam=1000, startModel=startModel)
print(model1)
print(inv.chi2(), inv.relrms(), pg.sum(inv.model[:nlay-1]))


# %%%
# To formulate the constraints, we need to set up a matrix for the left side
# and a vector for the right side of the equation. The layer thicknesses are
# the first values in the model vector, and setting 1 implements d1+d2+d3:
#      d1 d2 d3 r1 r2 r3 r4
# G = [1  1  1  0  0  0  0]
# c = [z2]
# The constraints G * m = c are set by setParameterConstraints with a
# Lagrangian parameter that determines how well the equation needs to be fit.
#

G = pg.Matrix(rows=1, cols=len(startModel))
for i in range(3):
    G.setVal(0, i, 1)

c = pg.Vector(1, pg.sum(synthk))
inv.setParameterConstraints(G, c, 100)
model2 = inv.run(rhoa, error, lam=1000, startModel=startModel)
print(model2)
print(inv.chi2(), inv.relrms(), pg.sum(inv.model[:nlay-1]))

# %%%
# All three models are plotted together and are equivalently fitting the data.
# Due to the additional constraints, the model is much closer to the synthetic.
#

fig, ax = pg.plt.subplots()
drawModel1D(ax, synthk, synres, plot="semilogx", label="synth")
drawModel1D(ax, model=model1, label="unconstrained")
drawModel1D(ax, model=model2, label="constrained")
ax.set_ylim(15, 0)
ax.grid(True)
ax.legend()

# %%%
# References
# ----------
# Wagner, F.M., Mollaret, C., Günther, T., Kemna, A., Hauck, A. (2019):
#     Quantitative imaging of water, ice, and air in permafrost systems through
#     petrophysical joint inversion of seismic refraction and electrical
#     resistivity data. Geophys. J. Int. 219, 1866-1875.
#     doi:10.1093/gji/ggz402.
