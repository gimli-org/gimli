#!/usr/bin/env python
# coding: utf-8
"""
Crosshole traveltime tomography
-------------------------------

Seismic and ground penetrating radar (GPR) methods are frequently applied to
image the shallow subsurface. While novel developments focus on inverting the
full waveform, ray-based approximations are still widely used in practice and
offer a computationally efficient alternative. Here we demonstrate the modeling
of traveltimes and their inversion for the underlying slowness distribution for
a  crosshole scenario.

We start by importing the necessary packages.
"""
# sphinx_gallery_thumbnail_number = 3

from itertools import product

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import pygimli as pg
import pygimli.meshtools as mt
from pygimli.physics.traveltime import Refraction

mpl.rcParams['image.cmap'] = 'inferno_r'


################################################################################
# Next, we build the crosshole acquisition geometry with two shallow boreholes.

# Acquisition parameters
bh_spacing = 20.0
bh_length = 25.0
sensor_spacing = 2.5

world = mt.createRectangle(start=[0, -(bh_length + 3)], end=[bh_spacing, 0.0],
                           marker=0)

depth = -np.arange(sensor_spacing, bh_length + sensor_spacing, sensor_spacing)

sensors = np.zeros((len(depth) * 2, 2))  # two boreholes
sensors[len(depth):, 0] = bh_spacing  # x
sensors[:, 1] = np.hstack([depth] * 2)  # y

################################################################################
# Traveltime calculations work on unstructured meshes and structured grids. We
# demonstrate this here by simulating the synthetic data on an unstructured mesh
# and inverting it on a simple structured grid.

# Create forward model and mesh
c0 = mt.createCircle(pos=(7.0, -10.0), radius=3, segments=25, marker=1)
c1 = mt.createCircle(pos=(12.0, -18.0), radius=4, segments=25, marker=2)
geom = mt.mergePLC([world, c0, c1])
for sen in sensors:
    geom.createNode(sen)

mesh_fwd = mt.createMesh(geom, quality=34, area=.25)
model = np.array([2000., 2300, 1700])[mesh_fwd.cellMarkers()]
pg.show(mesh_fwd, model, label="Velocity (m/s)", nLevs=3)

################################################################################
# Create inversion mesh
refinement = 0.25
x = np.arange(0, bh_spacing + refinement, sensor_spacing * refinement)
y = -np.arange(0.0, bh_length + 3, sensor_spacing * refinement)
mesh = pg.createMesh2D(x, y)

ax, _ = pg.show(mesh, hold=True)
ax.plot(sensors[:, 0], sensors[:, 1], "ro")

################################################################################
# Next, we create an empty DataContainer and fill it with sensor positions and
# all possible shot-recevier pairs for the two-borehole scenario using the
# prouct funtion in the itertools module (Python standard library).

numbers = np.arange(len(depth))
rays = list(product(numbers, numbers + len(numbers)))

# Empty container
scheme = pg.DataContainer()

# Add sensors
for sen in sensors:
    scheme.createSensor(sen)

# Add measurements
rays = np.array(rays)
scheme.resize(len(rays))
scheme.add("s", rays[:, 0])
scheme.add("g", rays[:, 1])
scheme.add("valid", np.ones(len(rays)))
scheme.registerSensorIndex("s")
scheme.registerSensorIndex("g")

################################################################################
# The forward simulation is performed with a few lines of code. We initialize an
# instance of the Refraction manager and call its `simulate` function with the
# mesh, the scheme and the slowness model (1 / velocity). We also add 0.1%
# relative and 10 microseconds of absolute noise.
#
# Secondary nodes allow for more accurate forward simulations. Check out the
# paper by `Giroux & Larouche (2013)
# <https://doi.org/10.1016/j.cageo.2012.12.005>`_ to learn more about it.

tt = Refraction()
mesh_fwd.createSecondaryNodes(5)
data = tt.simulate(mesh=mesh_fwd, scheme=scheme, slowness=1. / model,
                   noiseLevel=0.001, noiseAbs=1e-5)

################################################################################
# For the inversion we create a new instance of the Refraction manager to avoid
# confusion, since it is working on a different mesh.

ttinv = Refraction()
ttinv.setData(data)  # Set previously simulated data
ttinv.setMesh(mesh.createMeshWithSecondaryNodes(5))
invmodel = ttinv.invert(lam=3000, vtop=2000, vbottom=2000)
print("chi^2 =", ttinv.inv.getChi2())  # Look at the data fit

################################################################################
# Finally, we visualize the true model and the inversion result next to each
# other.
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 8), sharex=True, sharey=True)
ax1.set_title("True model")
ax2.set_title("Inversion result")

pg.show(mesh_fwd, model, ax=ax1, showMesh=True, label="Velocity (m/s)", nLevs=3)

for ax in (ax1, ax2):
    ax.plot(sensors[:, 0], sensors[:, 1], "wo")

ttinv.showResult(ax=ax2)
ttinv.showRayPaths(ax=ax2, color="0.8", alpha=0.3)
fig.tight_layout()

################################################################################
# Note how the rays are attracted by the high velocity anomaly while 
# circumventing the low velocity region.
