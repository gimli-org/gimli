# *-* coding: utf-8 *-*
# test some characteristics of ERTManager.simulate:
# 1) for complex conductivity models the response should not depend on the sign
#    of the K-factor

import numpy as np
import pygimli as pg
import pygimli.meshtools as mt
from pygimli.physics import ert

world = mt.createWorld(
    start=[-50, 0], end=[50, -50], layers=[-1, -5], worldMarker=True)

scheme = pg.DataContainerERT()
elec_positions = [[-2, 0, 0],
                  [-1, 0, 0],
                  [+0, 0, 0],
                  [+1, 0, 0],
                  [+2, 0, 0]]
for x, y, z in elec_positions:
    scheme.createSensor((x, y, z))

# Define number of measurements
scheme.resize(3)

# first two measurements have reversed geometric factor!
scheme.set('a', [0, 0, 3])
scheme.set('b', [1, 1, 2])
scheme.set('m', [3, 2, 1])
scheme.set('n', [2, 3, 0])

for pos in scheme.sensorPositions():
    world.createNode(pos)
    # world.createNode(pos + pg.RVector3(0, -0.1))

mesh = mt.createMesh(world, quality=34)

rhomap = [[1, 99.595 + 8.987j],
          [2, 99.595 + 8.987j],
          [3, 59.595 + 8.987j]]

# mgr = pg.physics.ERTManager()  # not necessary anymore
data = ert.simulate(mesh, res=rhomap, scheme=scheme, verbose=True)
rhoa = data.get('rhoa').array()
phia = data.get('phia').array()

# make sure all computed responses are equal, especially the first two, which
# only differ in the sign of their geometrical factor
np.testing.assert_allclose(rhoa, rhoa[0])
np.testing.assert_allclose(phia, phia[0])

# make sure the phases are positive (for positive input)
assert np.all(phia > 0)

# make sure rhoa is also positive
assert np.all(rhoa > 0)
