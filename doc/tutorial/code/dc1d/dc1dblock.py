#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import pygimli as pg
from pygimli.mplviewer import drawModel1D

nlay = 4
lam = 200.
errPerc = 3.

abmnr = pg.RMatrix()
pg.loadMatrixCol(abmnr, "sond1-100.ves")
ab2 = abmnr[0]
mn2 = abmnr[1]
rhoa = abmnr[2]

transRho = pg.RTransLogLU(1., 1000.)
transThk = pg.RTransLog()
transRhoa = pg.RTransLog()

f = pg.DC1dModelling(nlay, ab2, mn2)
f.region(0).setTransModel(transThk)
f.region(1).setTransModel(transRho)

paraDepth = max(ab2) / 3
f.region(0).setStartValue(max(ab2) / 3. / nlay / 2.)
f.region(1).setStartValue(pg.median(rhoa))

model = f.createStartVector()
model[nlay] *= 1.5

inv = pg.RInversion(rhoa, f, True)
inv.setModel(model)
inv.setTransData(transRhoa)
inv.setRelativeError(errPerc / 100.0)
inv.setLambda(lam)
inv.setMarquardtScheme(0.9)
model = inv.run()

fig, ax = plt.subplots(nrows=2)
ax[0].loglog(rhoa, ab2, 'rx-', inv.response(), ab2, 'b-')
ax[0].set_ylim(max(ab2), min(ab2))
ax[0].grid(which='both')
ax[0].set_xlabel(r"\rho_a in \Omegam")
ax[0].set_ylabel("AB/2 in m")
ax[0].set_legend(("measured", "fitted"), loc="upper left")
res, thk = model(nlay-1, nlay*2-1), model(0, nlay-1)
drawModel1D(ax[1], thk, res, xlabel=r'\rho in \Omega m', plot='semilogx')
drawModel1D(ax[1], [0.5, 3.5, 6.], [100., 500., 20., 1000.])
ax[1].grid(which='both')
plt.show()
