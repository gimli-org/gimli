"""Simple example of synthetic TDIP modelling and inversion."""
# %%
import numpy as np
import pygimli as pg
from pygimli.physics import ert
import pygimli.meshtools as mt
# %%
world = mt.createWorld(start=[-50, 0], end=[50, -50],
                       layers=[-1, -5], worldMarker=True)
scheme = ert.createData(
                    elecs=pg.utils.grange(start=-10, end=10, n=21),
                    schemeName='dd')
world += mt.createCircle(pos=[0, -3], radius=1, marker=4)
for pos in scheme.sensorPositions():
    _= world.createNode(pos)
    _= world.createNode(pos + [0.0, -0.1])
mesh = mt.createMesh(world, quality=34)
# %%
rhomap = [
   [1, 100. + 0j],
   [2, 50. + 0j],
   [3, 10.+ 0j],
   [4, 50.+ 1j],
]
dataFD = ert.simulate(mesh, res=rhomap, scheme=scheme, verbose=True)
dataFD.show("phia")
# %%
res = np.array([0, 100, 50, 10, 50.])
m = np.array([0, 0, 0, 0, 0.1])
mgr = ert.ERTIPManager()
dataTD = mgr.simulate(mesh=mesh, scheme=scheme, res=res, m=m)
# %%
dataTD.show("ip")
# %%
dataTD["err"] = 0.03
mgr = ert.ERTIPManager(dataTD)
mgr.invert(zWeight=0.2)
mgr.showResults()
