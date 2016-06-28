#!/usr/bin/env python
# encoding: utf-8

r"""
Meshing the Omega aka. BERT logo
================================

This is a fun example creating a logo for the BERT software. It illustrates the
possibility to hand over matplotlib path objects to the TriangleWrapper."""

import matplotlib.pyplot as plt
import matplotlib.textpath
import pygimli as pg

###############################################################################
# We start by generating a matplotlib path respresenting the :math:`\Omega`
# character.

logo_path = matplotlib.textpath.TextPath((0, 0), r'$\Omega$', size=1)
patch = matplotlib.patches.PathPatch(logo_path)

###############################################################################
# The vertices of the path are defined as mesh nodes and connected with edges.

nodes = patch.get_verts() * 50
nodes = pg.utils.unique_rows(nodes)  # remove duplicate nodes
poly = pg.Mesh(2)

for node in nodes:
    poly.createNode(node[0], node[1], 0.0)

for i in range(poly.nodeCount() - 1):
    poly.createEdge(poly.node(i), poly.node(i + 1))

poly.createEdge(poly.node(poly.nodeCount() - 1), poly.node(0))

###############################################################################
# We call the TriangleWrapper to generate the mesh and set the x values as the
# data for a color transition.

tri = pg.TriangleWrapper(poly)
mesh = pg.Mesh(2)
tri.setSwitches('-pzeAfa5q33')
tri.generate(mesh)

data = []
for cell in mesh.cells():
    data.append(cell.center().x())

###############################################################################
# Last, we create a BERT caption, visualize the mesh and fine-tune the figure.

fig, ax = plt.subplots(figsize=(4, 3))
ax.axis('off')
offset = -10
t = ax.text(1.5, offset, 'BERT', fontsize=50, fontweight='bold')
pg.show(mesh, data, axes=ax, cmap='RdBu', logScale=False, showLater=True)
pg.show(mesh, axes=ax)
ax.set_ylim(offset, mesh.ymax())

pg.wait()
