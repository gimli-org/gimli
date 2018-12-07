#!/usr/bin/env python
# encoding: utf-8

r"""
Region markers
==============
When constructing complex geometries out of basic geometric shapes (e.g.,
circle, rectangle, ...) we need to be careful with the region markers and their
positions. This example shows how to use markerPositions to properly set region
markers.
"""

###############################################################################
# imports
import pygimli as pg
import pygimli.meshtools as mt

###############################################################################
# In this first part we naively combine objects and assign markers to them,
# expecting two regions of concentric rings with markers 1 and 2. Note how the
# outer ring is assigned the marker 0 in the figure, although we specified
# marker=1 for the larger circle?
circle_outer = mt.createCircle(
    pos=[0.0, 0.0],
    radius=3.0,
    marker=1
)

circle_inner = mt.createCircle(
    pos=[0.0, 0.0],
    radius=1.0,
    # area=.3,
    boundaryMarker=0,
    marker=2
    )

plc = mt.mergePLC([circle_outer, circle_inner])

ax, cb = pg.show(plc, marker=True, savefig="plc_naive.pdf")

###############################################################################
# the solution to this problem are the region marker, which define the marker
# value of the region they are placed in. By default all region markers are
# assigned the position (0,0,0), thereby overwriting each other (see dots in
# figure below). If no region marker is present in a region, a marker value of
# 0 is assigned.

fig = ax.get_figure()
for nr, marker in enumerate(plc.regionMarker()):
    print(
        'Position marker number {}:'.format(nr + 1),
        marker.x(), marker.y(), marker.z()
    )
    ax.scatter(marker.x(), marker.y(), s=(2 - nr) * 20)
ax.set_title('marker positions')
fig.savefig('plc_naive_marker_positions.pdf', bbox_inches='tight')

###############################################################################
# Let fix this issue by assigning region marker positions that are not
# overwritten by other objects when the geometries are merged:
circle_outer = mt.createCircle(
    pos=[0.0, 0.0],
    radius=3.0,
    marker=1,
    markerPosition=[2.95, 0.0],
)

circle_inner = mt.createCircle(
    pos=[0.0, 0.0],
    radius=1.0,
    # area=.3,
    boundaryMarker=0,
    marker=2,
    markerPosition=[0.95, 0.0],
    )

plc = mt.mergePLC([circle_outer, circle_inner])

ax, cb = pg.show(plc, marker=True, savefig="plc_fixed.pdf")

fig = ax.get_figure()
for nr, marker in enumerate(plc.regionMarker()):
    print(
        'Position marker number {}:'.format(nr + 1),
        marker.x(), marker.y(), marker.z()
    )
    ax.scatter(marker.x(), marker.y(), s=(2 - nr) * 20)
ax.set_title('marker positions')
fig.savefig('plc_fixed_marker_positions.pdf', bbox_inches='tight')
