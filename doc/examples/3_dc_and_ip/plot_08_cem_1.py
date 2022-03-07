#!/usr/bin/env python
# encoding: utf-8
r"""
CEM 1: Using the Complete Electrode Model using PyGimli
-------------------------------------------------------
"""

###############################################################################
# .. warning::
#
#     Please remember that by using CEM electrodes in a mesh you are actively
#     changing the electrical characteristics of this mesh due to the spatially
#     distributed conducting electrode surface.
#     Therefore be careful with adding too many CEM electrodes in your mesh.
#
# .. warning::
#
#     Point electrodes can be used in conjunction with CEM electrodes and will
#     not influence the modeling results.
#     They can therefore be used as ideal potential probes.
#     HOWEVER, due to the specific implementation of free electrodes (i.e.,
#     electrodes that do not lie exactly on a mesh node AND are marked as nodes)
#     these free electrodes effectively are CEM electrodes in the sense that they
#     distribute current injection over the nearest mesh cell.
#     This can lead to potentially significant changes in modeling results.
#     To prevent this always add nodes with the corresponding marker for point
#     electrodes when CEM electrodes are also used. For example add an electrode
#     at x=1.15, y=0, z=0::
#
#         plc.createNode([1.15, 0, 0], marker=pg.core.MARKER_NODE_ELECTRODE)
#
#     Also make sure to add refining electrodes in the vicinity of any such point
#     electrode to make sure the mesh is fine enough::
#
#         plc.createNode([1.16, 0, 0])
#
#

###############################################################################
# * bug where calcK deletes contact impedances
# * Add section on extracting potentials on CEM electrodes, and thus how to
#   simulate 2-Point measurements
# * Plot potential distribution across extended potential electrodes
# * Investigate and document the effect of contact resistance and impedance;
#  show how to change using Python-interface
# * Add one series-circuit check

# Imports for this example
import numpy as np
import pygimli as pg
import pygimli.meshtools as mt
import pygimli.physics.ert as ert

import matplotlib.pylab as plt

###############################################################################
# A simple block with four CEM electrodes
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# We create a simple homogeneous 3D block, assign a resistivity and compute one
# four-point measurement and its reciprocal.
# Finally we provide the modeled measurement voltage at electrode M and N, as
# well as the apparent resistivity.
#
# .. figure:: ../../_static/cem/setup_only_electrodes.png
#     :align: center
#
#     ???

###############################################################################
# 1. Mesh generation
cube_width_x = 6
cube_height_y = 1
cube_length_z = 1

x_pot_electrode = 0.1
z_pot_electrode = 0.1

plc = mt.createCube(
    start=[0, 0, 0],
    end=[cube_width_x, cube_height_y, cube_length_z],
)

plc.createNode(
    [cube_width_x / 2, cube_height_y / 2, cube_length_z / 2],
    marker=-999
)
plc.createNode(
    [cube_width_x / 2, cube_height_y / 3, cube_length_z / 3],
    marker=-1000
)

# add corner nodes for CEM potential electrodes
offset = x_pot_electrode / 2
pot_electrode_edges = []
for nr, x_center in enumerate((cube_width_x / 3, cube_width_x / 3 * 2)):
    xmin = x_center - x_pot_electrode / 2 - offset
    xmax = x_center + x_pot_electrode / 2 + offset
    z_center = cube_length_z / 2
    zmin = z_center - z_pot_electrode / 2 - offset
    zmax = z_center + z_pot_electrode / 2 + offset
    # for x in (xmin, xmax):
    #     for z in (zmin, zmax):
    for x in np.linspace(xmin, xmax, 15):
        for z in np.linspace(zmin, zmax, 15):
            node = [x, 0, z]
            plc.createNode(node)
    pot_electrode_edges.append((
        xmin + offset,
        xmax - offset,
        zmin + offset,
        zmax - offset
    ))


mesh_coarse = mt.createMesh(
    plc,
    quality=32,
    # area=0.005,
)

# # additional refinements
mesh = mesh_coarse.createH2()
# mesh = mesh_coarse

# # Create a P2-optimized mesh (quadratic shape functions)
# mesh = mesh.createP2()

print(mesh)

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
# CEM electrodes are denoted by markers <= -10000
el_nr = 2
for (xmin, xmax, zmin, zmax) in pot_electrode_edges:
    for bound in mesh.boundaries():
        x = bound.center().x()
        y = bound.center().y()
        z = bound.center().z()
        if (x >= xmin and x <= xmax and y == 0 and z >= zmin and z <= zmax):
            bound.setMarker(-10000 - el_nr)
    el_nr += 1

# print(mesh)
# export for visualization using Paraview
mesh.exportVTK('mesh.vtk')
mesh.exportBoundaryVTU("boundaries.vtu")

###############################################################################
# 2. Create the measurement scheme
scheme = pg.DataContainerERT()

# electrode positions
pos = np.array((
    (0, 0, 0),
    (1, 0, 0),
    (2, 0, 0),
    (3, 0, 0),
))

# logical measurements
scheme.setSensorPositions(pos)
measurements = np.array((
    [0, 1, 2, 3],
    [2, 3, 0, 1],
))

for i, elec in enumerate("abmn"):
    scheme[elec] = measurements[:, i]

print('Computing geometric factors for later')
k = ert.createGeometricFactors(
    scheme, numerical=True, mesh=mesh, h2=False, p2=False
)

###############################################################################
# 3. Assign the resistivity model, here 10 k Ohm m

rhomap = [
    [1, 10000],
]
rho = pg.solver.parseArgToArray(rhomap, mesh.cellCount(), mesh)

###############################################################################
# Do the actual forward modeling
data = ert.simulate(
    mesh,
    res=rhomap,
    scheme=scheme,
    # no singularity removal for closed 3D geometries!
    sr=False,
    verbose=True,
    complex=False,
    # this makes sure we get the potentials. Turning this to False will yield
    # only a rhoa field in the returned data container
    calcOnly=True,
    # contactImpedances=[10, 10, 10, 10],
    # contactResistances=[100, 1000, 50, 10],
)

print(data)
print('--------------------')
print('Voltages M-N [V]:')
print(data['u'])

print('Transfer resistances [Ohm]')
print(data['u'] / data['i'])

print('Geometrical factors:')
print(k)

print('Apparent resistivities:')
# compute apparent resistivity by hand
data['rhoa'] = data['u'] / data['i'] * k
print(data['rhoa'])

print('Expected [Ohm m]:', rhomap[0][1])

###############################################################################
# Lets examine the potentials across the potential electrodes for very small
# and very high contact resistances

node_potentials_lowc = ert.simulate(
        mesh,
        res=rhomap,
        scheme=scheme,
        sr=False,
        verbose=True,
        complex=False,
        # this flag will return the node potentials
        returnFields=True,
        contactImpedances=[1e-6, 1e-6, 1e-6, 1e-6],
)

node_potentials_highc = ert.simulate(
        mesh,
        res=rhomap,
        scheme=scheme,
        sr=False,
        verbose=True,
        complex=False,
        # this flag will return the node potentials
        returnFields=True,
        contactImpedances=[1e6, 1e6, 1e6, 1e6],
)

# extract node potentials for the fourth CEM electrode

# potential electrodes have ids -10002 and -10003
r = mesh.findBoundaryByMarker(-10003)
nodes = np.unique(
    np.vstack([np.vstack([item.id() for item in q.nodes()]) for q in r])
)
node_positions = mesh.positions(nodes).array()

node_pots = node_potentials_lowc[2][nodes].array()
potentials_lowc = np.hstack((node_positions, node_pots[:, np.newaxis]))

node_pots = node_potentials_highc[2][nodes].array()
potentials_highc = np.hstack((node_positions, node_pots[:, np.newaxis]))

# cem_nodes = np.where(np.array(mesh.boundaryMarkers()) == -10002)[0]
# # electrode_potentials = node_potentials[
# positions = mesh.positions(cem_nodes).array()
# pots = node_potentials[2][cem_nodes].array()

fig, axes = plt.subplots(1, 2)
ax = axes[0]
ax.set_title('Low contact resistances')
sc = ax.scatter(
    potentials_lowc[:, 0],
    potentials_lowc[:, 2],
    c=potentials_lowc[:, 3],
    s=500,
    cmap="jet",
    vmin=29,
    vmax=35,
)
cb = fig.colorbar(sc, ax=ax)
cb.set_label('Potential [V]')
ax.set_ylabel('z [m]')
ax.set_xlabel('x [m]')
ax = axes[1]
ax.set_title('High contact resistances')
sc = ax.scatter(
    potentials_highc[:, 0],
    potentials_highc[:, 2],
    c=potentials_highc[:, 3],
    s=500,
    cmap="jet",
    vmin=29,
    vmax=35,
)
cb = fig.colorbar(sc, ax=ax)
ax.set_xlabel('x [m]')
cb.set_label('Potential [V]')
fig.show()

###############################################################################
# The effect of the contact impedances is clearly visible: the quasi-short
# circuit for low contact impedances leads to an equal potential across the
# electrode surface, whereas a high contact impedance will lead to a potential
# drop across the electrode surface.
#
# .. note::
#
#   The node potentials can also be used to simulate 2-Point measurements. Just
#   extract note potentials for the current electrodes and average those
#   potentials to get a simlulated 2P voltage measurement for an injection
#   current of 1 A.

###############################################################################
# Validation 1: Two blocks in series
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# TODO: TEXT

###############################################################################
# Mesh generation
cube_width_x = 3
cube_height_y = 1
cube_length_z = 1

# size of potential electrodes
x_pot_electrode = 0.01
z_pot_electrode = 0.01

# we construct subvolumes: 4x4 in x-y plane, and 4 in z
Nx = 3
Ny = 1
Nz = 1

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

plc.createNode(
    [cube_width_x / 2, cube_height_y / 2, cube_length_z / 2],
    marker=-999
)
plc.createNode(
    [cube_width_x / 2, cube_height_y / 3, cube_length_z / 3],
    marker=-1000
)

# add corner nodes for CEM potential electrodes
offset = x_pot_electrode / 2
pot_electrode_edges = []
for nr, x_center in enumerate((cube_width_x / 3, cube_width_x / 3 * 2)):
    xmin = x_center - x_pot_electrode / 2 - offset
    xmax = x_center + x_pot_electrode / 2 + offset
    z_center = cube_length_z / 2
    zmin = z_center - z_pot_electrode / 2 - offset
    zmax = z_center + z_pot_electrode / 2 + offset
    # for x in (xmin, xmax):
    #     for z in (zmin, zmax):
    for x in np.linspace(xmin, xmax, 15):
        for z in np.linspace(zmin, zmax, 15):
            node = [x, 0, z]
            plc.createNode(node)
    pot_electrode_edges.append((
        xmin + offset,
        xmax - offset,
        zmin + offset,
        zmax - offset
    ))


mesh_coarse = mt.createMesh(
    plc,
    quality=32,
    # area=0.005,
)

# # additional refinements
mesh = mesh_coarse.createH2()
# mesh = mesh_coarse

# # Create a P2-optimized mesh (quadratic shape functions)
# mesh = mesh.createP2()

print(mesh)

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
# CEM electrodes are denoted by markers <= -10000
el_nr = 2
for (xmin, xmax, zmin, zmax) in pot_electrode_edges:
    for bound in mesh.boundaries():
        x = bound.center().x()
        y = bound.center().y()
        z = bound.center().z()
        if (x >= xmin and x <= xmax and y == 0 and z >= zmin and z <= zmax):
            bound.setMarker(-10000 - el_nr)
    el_nr += 1
mesh.exportVTK('mesh_series.vtk')
mesh.exportBoundaryVTU("boundaries_series.vtu")

###############################################################################
# 2. Create the measurement scheme
scheme = pg.DataContainerERT()

# electrode positions
pos = np.array((
    (0, 0, 0),
    (1, 0, 0),
    (2, 0, 0),
    (3, 0, 0),
))

# logical measurements
scheme.setSensorPositions(pos)
measurements = np.array((
    [0, 1, 2, 3],
    [2, 3, 0, 1],
))

for i, elec in enumerate("abmn"):
    scheme[elec] = measurements[:, i]

scheme.set('k', [1 for x in range(scheme.size())])
# print('Computing geometric factors for later')
# k = ert.createGeometricFactors(
#     scheme, numerical=True, mesh=mesh, h2=False, p2=False
# )

###############################################################################
# 3. Assign the resistivity model, here 10 k Ohm m

rhomap = [
    [1, 1000],
    [2, 2000],
    [3, 1000],
]
rho = pg.solver.parseArgToArray(rhomap, mesh.cellCount(), mesh)

###############################################################################
# Do the actual forward modeling
node_potentials = ert.simulate(
    mesh,
    res=rhomap,
    scheme=scheme,
    # no singularity removal for closed 3D geometries!
    sr=False,
    verbose=True,
    # calcOnly=True,
    returnFields=True,
    contactImpedances=[0.1, 0.1, 0.1, 0.1],
    # contactResistances=[100, 1000, 50, 10],
)

# print(node_potentials)
# print('--------------------')
# print('Transfer resistances [Ohm]')
# print(node_potentials['u'] / node_potentials['i'])

# print('Apparent resistivities:')
# # compute apparent resistivity by hand
# data['rhoa'] = data['u'] / data['i'] * k
# print(data['rhoa'])

# superpositions of electrode 0 and electrode 1
first_injection = node_potentials[0] - node_potentials[1]


def get_cem_voltage(marker):
    # return the averaged electrode potential
    r = mesh.findBoundaryByMarker(marker)
    nodes = np.unique(
        np.vstack([np.vstack([item.id() for item in q.nodes()]) for q in r])
    )
    # node_positions = mesh.positions(nodes).array()
    node_pots = first_injection[nodes].array()
    # potentials = np.hstack((node_positions, node_pots[:, np.newaxis]))
    return np.mean(node_pots)


rho1 = rhomap[0][1]
rho2 = rhomap[1][1]
rho3 = rhomap[2][1]

le = cube_width_x / 3
A = cube_height_y * cube_length_z

Z1 = rho1 * A / le
Z2 = rho2 * A / le
Z3 = rho3 * A / le
Zg = Z1 + Z2 + Z3

voltage_2p = get_cem_voltage(-10000) - get_cem_voltage(-10001)

print('Transfer resistance [Ohm]:', voltage_2p / 1.0)
print('Expected transfer resistance [Ohm]:', Zg)

###############################################################################
# Validation 2: Two blocks in parallel
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# TODO: TEXT


###############################################################################
# Mesh generation
cube_width_x = 1
cube_height_y = 1
cube_length_z = 1

# size of potential electrodes
x_pot_electrode = 0.01
z_pot_electrode = 0.01

# we construct subvolumes: 4x4 in x-y plane, and 4 in z
Nx = 1
Ny = 2
Nz = 1

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

plc.createNode(
    [cube_width_x / 2, cube_height_y / 2, cube_length_z / 2],
    marker=-999
)
plc.createNode(
    [cube_width_x / 2, cube_height_y / 3, cube_length_z / 3],
    marker=-1000
)

# add corner nodes for CEM potential electrodes
offset = x_pot_electrode / 2
pot_electrode_edges = []
for nr, x_center in enumerate((cube_width_x / 3, cube_width_x / 3 * 2)):
    xmin = x_center - x_pot_electrode / 2 - offset
    xmax = x_center + x_pot_electrode / 2 + offset
    z_center = cube_length_z / 2
    zmin = z_center - z_pot_electrode / 2 - offset
    zmax = z_center + z_pot_electrode / 2 + offset
    # for x in (xmin, xmax):
    #     for z in (zmin, zmax):
    for x in np.linspace(xmin, xmax, 15):
        for z in np.linspace(zmin, zmax, 15):
            node = [x, 0, z]
            plc.createNode(node)
    pot_electrode_edges.append((
        xmin + offset,
        xmax - offset,
        zmin + offset,
        zmax - offset
    ))


mesh_coarse = mt.createMesh(
    plc,
    quality=32,
    # area=0.005,
)

# # additional refinements
mesh = mesh_coarse.createH2()
# mesh = mesh_coarse

# # Create a P2-optimized mesh (quadratic shape functions)
# mesh = mesh.createP2()

print(mesh)

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
# CEM electrodes are denoted by markers <= -10000
el_nr = 2
for (xmin, xmax, zmin, zmax) in pot_electrode_edges:
    for bound in mesh.boundaries():
        x = bound.center().x()
        y = bound.center().y()
        z = bound.center().z()
        if (x >= xmin and x <= xmax and y == 0 and z >= zmin and z <= zmax):
            bound.setMarker(-10000 - el_nr)
    el_nr += 1
mesh.exportVTK('mesh_series.vtk')
mesh.exportBoundaryVTU("boundaries_series.vtu")

###############################################################################
# 2. Create the measurement scheme
scheme = pg.DataContainerERT()

# electrode positions
pos = np.array((
    (0, 0, 0),
    (1, 0, 0),
    (2, 0, 0),
    (3, 0, 0),
))

# logical measurements
scheme.setSensorPositions(pos)
measurements = np.array((
    [0, 1, 2, 3],
    [2, 3, 0, 1],
))

for i, elec in enumerate("abmn"):
    scheme[elec] = measurements[:, i]

scheme.set('k', [1 for x in range(scheme.size())])
# print('Computing geometric factors for later')
# k = ert.createGeometricFactors(
#     scheme, numerical=True, mesh=mesh, h2=False, p2=False
# )

###############################################################################
# 3. Assign the resistivity model, here 10 k Ohm m

rhomap = [
    [1, 100],
    [2, 5000],
]
rho = pg.solver.parseArgToArray(rhomap, mesh.cellCount(), mesh)

###############################################################################
# Do the actual forward modeling
node_potentials = ert.simulate(
    mesh,
    res=rhomap,
    scheme=scheme,
    # no singularity removal for closed 3D geometries!
    sr=False,
    verbose=True,
    # calcOnly=True,
    returnFields=True,
    contactImpedances=[0.1, 0.1, 0.1, 0.1],
    # contactResistances=[100, 1000, 50, 10],
)

# print(node_potentials)
# print('--------------------')
# print('Transfer resistances [Ohm]')
# print(node_potentials['u'] / node_potentials['i'])

# print('Apparent resistivities:')
# # compute apparent resistivity by hand
# data['rhoa'] = data['u'] / data['i'] * k
# print(data['rhoa'])

# superpositions of electrode 0 and electrode 1
first_injection = node_potentials[0] - node_potentials[1]


def get_cem_voltage(marker):
    # return the averaged electrode potential
    r = mesh.findBoundaryByMarker(marker)
    nodes = np.unique(
        np.vstack([np.vstack([item.id() for item in q.nodes()]) for q in r])
    )
    # node_positions = mesh.positions(nodes).array()
    node_pots = first_injection[nodes].array()
    # potentials = np.hstack((node_positions, node_pots[:, np.newaxis]))
    return np.mean(node_pots)


rho1 = rhomap[0][1]
rho2 = rhomap[1][1]

le = cube_width_x
A = cube_height_y * cube_length_z / Ny
print('A', A)
print('le', le)

Z1 = rho1 * le / A
Z2 = rho2 * le / A
print('Z1, Z2', Z1, Z2)
Zg = 1 / (1 / Z1 + 1 / Z2)

voltage_2p = get_cem_voltage(-10000) - get_cem_voltage(-10001)

print('Transfer resistance [Ohm]:', voltage_2p / 1.0)
print('Expected transfer resistance [Ohm]:', Zg)
