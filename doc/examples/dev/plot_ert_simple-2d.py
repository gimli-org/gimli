#!/usr/bin/env python
# encoding: UTF-8

"""
    Super simple but complete 2D ERT without topography
"""

import pygimli as pg
import pybert as pb

data = pb.DataContainerERT('gallery.data')
print(data)

#mesh = pg.meshtools.createParaMesh2dGrid(data.sensorPositions(),
                                         #paraDZ=0.5)
mesh = pg.meshtools.createParaMesh(data.sensorPositions(), verbose=1, paraMaxCellSize=1,
                                   quality=34, smooth=[1,4], paraDepth=10, paraDX=0.5)

print(mesh)

#fop = pb.DCSRMultiElectrodeModelling(mesh, data)
fop = pg.physics.ert.ERTModelling(mesh, data)
fop.regionManager().region(1).setBackground(True)
fop.createRefinedForwardMesh(refine=True, pRefine=False)

inv = pg.Inversion(data("rhoa"), fop, verbose=True, dosave=True)

datTrans = pg.trans.TransLog()
modTrans = pg.trans.TransLog()

inv.setMaxIter(10)
inv.setTransData(datTrans)
inv.setTransModel(modTrans)
inv.setError(data('err'))
inv.setModel(pg.Vector(fop.regionManager().parameterCount(),
                        pg.math.median(data('rhoa'))))
inv.setLambda(5)

model = inv.run()

modelMesh = fop.regionManager().paraDomain()
a, cbar = pg.show(modelMesh, model)

##C = fop.contraintsMatrix()
#S = fop.jacobian()

#print(S)
#modelMesh = fop.regionManager().paraDomain()

#fig = pg.plt.figure()
#ax = [fig.add_subplot(2,2,i) for i in range(1,5)]
#a, cbar = pg.show(modelMesh, model, axes=ax[0])
#a, cbar = pg.show(modelMesh, model, tri=1, interpolate=1, axes=ax[1])
#a, cbar = pg.show(modelMesh, model, tri=1, interpolate=0, shading='gouraud', axes=ax[2])

#sens = pb.prepExportSensitivityData(modelMesh, S[110])
#a, cbar = pg.show(modelMesh, sens, axes=ax[3], cMap='b2r')

#for i in range(4):
    #pg.show(modelMesh, data.sensorPositions(), axes=ax[i])
    #pg.show(modelMesh, axes=ax[i])

#potmatrix = fop.solution()


#print(potmatrix)
#print(potmatrix[0])

#u = potmatrix[0] - potmatrix[1]

#uc = pg.Vector(mesh.cellCount())
#for c in mesh.cells():
    #uc[c.id()] = c.pot(c.center(), u)

#fig = pg.plt.figure()
#ax = [fig.add_subplot(3,2,i) for i in range(1,3*2+1)]

#pg.show(fop.mesh(), pg.logTransDropTol(u, 1e-5), tri=0, axes=ax[0])
#pg.show(fop.mesh(), pg.logTransDropTol(u, 1e-5), tri=1, axes=ax[1])
#pg.show(mesh, pg.logTransDropTol(uc, 1e-5), tri=0, axes=ax[2])
#pg.show(mesh, pg.logTransDropTol(uc, 1e-5), tri=1, interpolate=0, axes=ax[3])
#pg.show(mesh, pg.logTransDropTol(uc, 1e-5), tri=1, interpolate=1, axes=ax[4],
        #nLevs=256, omitLines=1, )
#pg.show(mesh, pg.logTransDropTol(uc, 1e-5), tri=1, axes=ax[5],
        #shading='gouraud', )

pg.wait()
