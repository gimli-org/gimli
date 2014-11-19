#!/usr/bin/env python
# encoding: UTF-8

"""
    Super simple but complete 2D ERT without topography
"""

import pygimli as pg
import pybert as pb

data = pb.DataContainerERT('gallery.dat')
print(data)

mesh = pg.meshtools.createParaMesh2dGrid(data.sensorPositions(),
                                         paraDZ=0.5)
print(mesh)

fop = pb.DCSRMultiElectrodeModelling(mesh, data)
fop.regionManager().region(1).setBackground(True)
fop.createRefinedForwardMesh(refine=True, pRefine=False)

inv = pg.RInversion(data("rhoa"), fop, verbose=True, dosave=True)

datTrans = pg.RTransLog()
modTrans = pg.RTransLog()

inv.setTransData(datTrans)
inv.setTransModel(modTrans)
inv.setError(data('err'))
inv.setModel(pg.RVector(fop.regionManager().parameterCount(),
                        pg.median(data('rhoa'))))
inv.setLambda(5)

model = inv.run()

modelMesh = fop.regionManager().paraDomain()

ax, cbar = pg.show(modelMesh, model, showLater=1)
pg.show(modelMesh, data.sensorPositions(), axes=ax, showLater=1)
pg.show(modelMesh, axes=ax)
