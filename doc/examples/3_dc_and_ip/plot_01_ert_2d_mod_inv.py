# -*- coding: utf-8 -*-
"""
2D ERT modeling and inversion
-----------------------------
"""

import numpy as np
import pygimli as pg
import pygimli.meshtools as mt
import pygimli.physics.ert as ert

###############################################################################
# Create geometry definition for the modelling domain.
# worldMarker=True indicates the default boundary conditions for the ERT
world = mt.createWorld(start=[-50, 0], end=[50, -50], layers=[-1, -5],
                       worldMarker=True)

###############################################################################
# Create some heterogeneous circular
block = mt.createCircle(pos=[0, -3.], radius=1, marker=4, boundaryMarker=10,
                        area=0.1)

###############################################################################
# Merge geometry definition into a Piecewise Linear Complex (PLC)
geom = mt.mergePLC([world, block])

###############################################################################
# Optional: show the geometry
pg.show(geom)

###############################################################################
# Create a Dipole Dipole ('dd') measuring scheme with 21 electrodes.
scheme = ert.createERTData(elecs=np.linspace(start=-10, stop=10, num=21),
                           schemeName='dd')

###############################################################################
# Put all electrodes (aka. sensors positions) into the PLC to enforce mesh
# refinement. Due to experience, its convenient to add further refinement
# nodes in a distance of 10% of electrode spacing, to achieve sufficient
# numerical accuracy.
for p in scheme.sensors():
    geom.createNode(p)
    geom.createNode(p - [0, 0.1])

# Create a mesh for the finite element modelling with appropriate mesh quality.
mesh = mt.createMesh(geom, quality=34)

# Create a map to set resistivity values in the appropriate regions
# [[regionNumber, resistivity], [regionNumber, resistivity], [...]
rhomap = [[1, 100.],
          [2, 50.],
          [3, 10.],
          [4, 100.]]

# Optional: take a look at the mesh
pg.show(mesh, data=rhomap, label=pg.unit('res'), showMesh=True)

###############################################################################
# Perform the modeling with the mesh and the measuring scheme itself
# and return a data container with apparent resistivity values,
# geometric factors and estimated data errors specified by the noise setting.
# The noise is also added to the data.
data = ert.simulate(mesh, res=rhomap, scheme=scheme, noiseLevel=1,
                    noiseAbs=1e-6)

pg.info('Simulated data', data)
pg.info('The data contains:', data.dataMap().keys())

pg.info('Simulated rhoa (min/max)', min(data['rhoa']), max(data['rhoa']))
pg.info('Selected data noise %(min/max)', min(data['err'])*100, max(data['err'])*100)

###############################################################################
# Optional: you can filter all values and tokens in the data container.
# Its possible that there are some negative data values due to noise and
# huge geometric factors. So we need to remove them
data.remove(data['rhoa'] < 0)
pg.info('Filtered rhoa (min/max)', min(data['rhoa']), max(data['rhoa']))

# Optional: you can save the data for further use
data.save('simple.dat')

# Optional: and you can take a look at the data
ert.show(data)

###############################################################################
# Initialize the ERTManager. You can initialize the manager with a data file 
# itself or a filename.
mgr = ert.ERTManager('simple.dat')

###############################################################################
# Run the inversion with the preset data. The Inversion mesh will be created 
# with default settings.
mgr.invert(lam=20, verbose=True)

###############################################################################
# Let the ERTManger show you the model and fitting results of the last
# successful run. Shows data, model response, and model.
#
mgr.showResultAndFit()

###############################################################################
# You can also provide your own mesh (e.g., a structured grid if you like them)
# 
inversionDomain = pg.createGrid(x=np.linspace(start=-12, stop=12, num=33),
                                y=-pg.cat([0], pg.utils.grange(0.5, 8, n=8)),
                                marker=2)
###############################################################################
# The inversion domain for ERT problems needs a boundary which represents the 
# far regions in the subsurface of the halfspace.
# Give a cell marker lower than the marker for the inversion region, the lowest 
# cell marker in the mesh will be the inversion boundary region.
#
mesh = pg.meshtools.appendTriangleBoundary(inversionDomain, marker=1,
                                           xbound=50, ybound=50)

###############################################################################
# The Inversion can be called with data and mesh as argument as well
#
model = mgr.invert(data, mesh=mesh, lam=20, verbose=True)

###############################################################################
# You can of course access to mesh and model and plot them for your own.
# Note, that the cells of the parametric domain of your mesh might be in 
# a different order than the values in the model array. 
# The manager can help to permutate them into the right order.
#
modelPD = mgr.paraModel(model)
pg.show(mgr.paraDomain, modelPD, label='Model', cMap='Spectral_r',
        logScale=True, cMin=10, cMax=100)

pg.info('Inversion stopped with chi² = {0:.3}'.format(mgr.fw.chi2()))
