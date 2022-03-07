#!/usr/bin/env python
r"""
CEM 2: Complex CEM Modeling using PyGimli
-----------------------------------------

TODO
====

* [ ] check output of the core, remove debug messages

Test cases
==========

* Series-circuit

  markers 1-32 == 100 * exp(1j * -50 / 1000)
  markers 33-64 == 50 * exp(1j * -30 / 1000)

  Zg = Z1 + Z2

  Z1 = rho1 * l / A
  Z2 = rho2 * l / A
  l = 2
  A = 1 x 1 = 1

"""
import numpy as np
import pygimli as pg
import pygimli.meshtools as mt
import pygimli.physics.ert as ert
###############################################################################

###############################################################################
# First we construct a mesh consisting of 64 regions, arranged in a 4x4x4
# geometry. The overall mesh shall have dimensions of
# * x: 1 meter
# * y: 1 meter
# * z: 4 meter
# Each region gets a unique marker ranging from 1 to 64, with the first 16
# lying in the x-y plane.
# We later use this mesh to test various complex conductivity parameterizations
# against analytical solutions.
#
cube_width_x = 4
cube_height_y = 1
cube_length_z = 1

# size of potential electrodes
x_pot_electrode = 0.1
z_pot_electrode = 0.1

# we construct subvolumes: 4x4 in x-y plane, and 4 in z
Nx = 4
Ny = 4
Nz = 4

increment_x = cube_width_x / Nx
increment_y = cube_height_y / Ny
increment_z = cube_length_z / Nz

cube_list = []
marker = 1
for i_x in range(Nx):
    for i_y in range(Ny):
        for i_z in range(Nz):
            start = [
                0 + i_x * increment_x,
                0 + i_y * increment_y,
                0 + i_z * increment_z,
            ]
            end = [
                0 + (i_x + 1) * increment_x,
                0 + (i_y + 1) * increment_y,
                0 + (i_z + 1) * increment_z,
            ]
            cube = mt.createCube(
                start=start,
                end=end,
                marker=marker,
            )
            marker += 1
            cube_list.append(cube)

plc = mt.mergePLC(cube_list)

# for closed geometries we need reference and calibration nodes, far away from
# any electrode
plc.createNode(
    [cube_width_x / 2, cube_height_y / 2, cube_length_z / 2],
    marker=pg.core.MARKER_NODE_REFERENCEELECTRODE,
)

plc.createNode(
    [cube_width_x / 2, cube_height_y / 3, cube_length_z / 3],
    marker=pg.core.MARKER_NODE_CALIBRATION
)

# add corner nodes for CEM potential electrodes
# this is not strictly necessary, but provides us with a nice regular node
# pattern across the CEM electrodes
pot_electrode_edges = []
for nr, x_center in enumerate((cube_width_x / 3, cube_width_x / 3 * 2)):
    xmin = x_center - x_pot_electrode / 2
    xmax = x_center + x_pot_electrode / 2
    z_center = cube_length_z / 2
    zmin = z_center - z_pot_electrode / 2
    zmax = z_center + z_pot_electrode / 2
    # for x in (xmin, xmax):
    #     for z in (zmin, zmax):
    for x in np.linspace(xmin, xmax, 5):
        for z in np.linspace(zmin, zmax, 5):
            node = [x, 0, z]
            plc.createNode(node)
    pot_electrode_edges.append((xmin, xmax, zmin, zmax))


mesh_coarse = mt.createMesh(
    plc,
    # quality=32,
)

# # additional refinements
mesh = mesh_coarse.createH2()

# # Create a P2-optimized mesh (quadratic shape functions)
# mesh = mesh.createP2()

# print(mesh)

# Set CEM boundary markers
# electrode at x = 0
for bound in mesh.boundaries():
    if np.isclose(bound.center().x(), 0):
        bound.setMarker(-10000)

# electrode at x = cube_width_x
for bound in mesh.boundaries():
    if np.isclose(bound.center().x(), cube_width_x):
        bound.setMarker(-10001)
# potential electrodes
el_nr = 2
for (xmin, xmax, zmin, zmax) in pot_electrode_edges:
    for bound in mesh.boundaries():
        x = bound.center().x()
        y = bound.center().y()
        z = bound.center().z()
        if (x >= xmin and x <= xmax and y == 0 and z >= zmin and z <= zmax):
            bound.setMarker(-10000 - el_nr)
    el_nr += 1

mesh.exportVTK('mesh.vtk')
mesh.exportBoundaryVTU("boundaries.vtu")

###############################################################################
# Create the measurement scheme
scheme = pg.DataContainerERT()

# we need to define electrodes in the scheme, even if some of them are CEM
# electrodes
# TODO: Explain the assignment scheme of PG here
pos = np.array((
    (0, 0, 0),
    (1, 0, 0),
    (2, 0, 0),
    (3, 0, 0),
))

scheme.setSensorPositions(pos)

# These are the logical four-point measurements, here a classical Wenner
# configuration such as is common for laboratory measurements. We also add the
# reciprocal measurement
measurements = np.array((
    [0, 1, 2, 3],
    [2, 3, 0, 1],
))

# register the measurements in the scheme object
for i, elec in enumerate("abmn"):
    scheme[elec] = measurements[:, i]

# in case we want to compute transfer resistances instead of apparent
# resistivities, set all geometric factors to 1
# scheme.set('k', [1 for x in range(scheme.size())])

###############################################################################
# Now we parameterize the model

rhomap = []
for i in range(1, 65):
    rhomap.append(
        [i, 100 - 1j]
    )
    # [1, 100 - 1j],
rho = pg.solver.parseArgToArray(rhomap, mesh.cellCount(), mesh)

# # for nr, c in enumerate(mesh.cells()):
# #     c.setMarker(nr)

###############################################################################
data = ert.simulate(
    mesh,
    res=rhomap,
    scheme=scheme,
    sr=False,
    verbose=True,
    complex=True,
    # noiseAbs=0.0,
    # noiseLevel=0.0,
    calcOnly=False,
)
print('Numerical:', data['k'])
# print('Theoretical:', compute_theoretical_k())
print(data['rhoa'])
print('phi', data['phia'])
# print('Expected:', compute_voltage())
