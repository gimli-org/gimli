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

inv.setMaxIter(1)
inv.setTransData(datTrans)
inv.setTransModel(modTrans)
inv.setError(data('err'))
inv.setModel(pg.RVector(fop.regionManager().parameterCount(),
                        pg.median(data('rhoa'))))
inv.setLambda(5)

model = inv.run()

# C = fop.contraintsMatrix()
# S = fop.jacobian()

modelMesh = fop.regionManager().paraDomain()
pg.showLater(1)

ax, cbar = pg.show(modelMesh, model)

pg.show(modelMesh, data.sensorPositions(), axes=ax, showLater=1)
pg.show(modelMesh, axes=ax)


ax, cbar = pg.show(modelMesh, model, tri=1, interpolate=1)
ax, cbar = pg.show(modelMesh, model, tri=1, interpolate=0, shading='gouraud')

potmatrix = fop.solution()
print(potmatrix)
print(potmatrix[0])

u = potmatrix[0] - potmatrix[1]
ax, cbar = pg.show(fop.mesh(), pg.logTransDropTol(u, 1e-5))

uc = pg.RVector(mesh.cellCount())
for c in mesh.cells():
    uc[c.id()] = c.pot(c.center(), u)

ax, cbar = pg.show(mesh, pg.logTransDropTol(uc, 1e-5))
ax, cbar = pg.show(mesh, pg.logTransDropTol(uc, 1e-5), tri=1, interpolate=0)
ax, cbar = pg.show(mesh, pg.logTransDropTol(uc, 1e-5), tri=1, interpolate=1,
                   nLevs=256, omitLines=1)
ax, cbar = pg.show(mesh, pg.logTransDropTol(uc, 1e-5), tri=1,
                   shading='gouraud')

pg.showNow()
