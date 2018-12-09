#!/usr/bin/env python
# encoding: utf-8

r"""
Region markers
==============

This example highlights the usage of regionMarkers, as well as some properties
of regions if the :py:module:`pygimli.meshtools` module is used to create
complex compound meshes.

When constructing complex geometries out of basic geometric shapes (e.g.,
circle, rectangle, ...) we need to be careful with the region markers and their
positions. This example shows how to use markerPositions to properly set region
markers.
"""

###############################################################################
# imports
import pygimli as pg
import pygimli.meshtools as mt

"""
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

ax, cb = pg.show(plc, markers=True)

###############################################################################
# The solution to this problem are the region marker, which define the marker
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
fig.show()

###############################################################################
# Let us fix this issue by assigning region marker positions that are not
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
    marker=2,
    markerPosition=[0.95, 0.0],
)

plc = mt.mergePLC([circle_outer, circle_inner])

ax, cb = pg.show(plc, markers=True)

fig = ax.get_figure()
for nr, marker in enumerate(plc.regionMarker()):
    print(
        'Position marker number {}:'.format(nr + 1),
        marker.x(), marker.y(), marker.z()
    )
    ax.scatter(marker.x(), marker.y(), s=(2 - nr) * 20)
ax.set_title('marker positions')
fig.show()
"""
###############################################################################
polygon1 = mt.createPolygon(
    [[0.0, 0.0], [1.0, 0.0], [1.0, -1.0], [0.0, -1]],
    isClosed=True,
    marker=1,
)

polygon2 = mt.createPolygon(
    [[0.25, -0.25], [0.75, -0.25], [0.75, -0.75], [0.25, -0.75]],
    isClosed=True,
    marker=2,
)

plc = mt.mergePLC([polygon1, polygon2])

ax, cb = pg.show(plc, markers=True)

fig = ax.get_figure()
for nr, marker in enumerate(plc.regionMarker()):
    print(
        'Position marker number {}:'.format(nr + 1),
        marker.x(), marker.y(), marker.z()
    )
    ax.scatter(marker.x(), marker.y(), s=(4 - nr) * 20, color='k')
ax.set_title('marker positions')
fig.savefig('polygons_test_abs.pdf')

###############################################################################
polygon1 = mt.createPolygon(
    [[0.0, 0.0], [1.0, 0.0], [1.0, -1.0], [0.0, -1]],
    isClosed=True,
    marker=1,
    markerPosition=[0.9, -0.9],
)

polygon2 = mt.createPolygon(
    [[0.25, -0.25], [0.75, -0.25], [0.75, -0.75], [0.25, -0.75]],
    isClosed=True,
    marker=2,
)

plc = mt.mergePLC([polygon1, polygon2])

ax, cb = pg.show(plc, markers=True)

fig = ax.get_figure()
for nr, marker in enumerate(plc.regionMarker()):
    print(
        'Position marker number {}:'.format(nr + 1),
        marker.x(), marker.y(), marker.z()
    )
    ax.scatter(marker.x(), marker.y(), s=(4 - nr) * 20, color='k')
ax.set_title('marker positions')
fig.savefig('polygons_test_abs_fixed.pdf')
###############################################################################
polygon1 = mt.createPolygon(
    [[0.0, 0.0], [1.0, 0.0], [1.0, -1.0], [0.0, -1]],
    isClosed=True,
    marker=1,
)

polygon2 = mt.createPolygon(
    [[1.0, -1.0], [2.0, -1.0], [2.0, -2.0], [1.0, -2.0]],
    isClosed=True,
    marker=2,
)

polygon3 = mt.createPolygon(
    [[0.6, -0.5], [1.4, -0.5], [1.4, -1.5], [0.6, -1.5]],
    isClosed=True,
    marker=3,
)

plc = mt.mergePLC([polygon1, polygon2, polygon3])

ax, cb = pg.show(plc, markers=True)

fig = ax.get_figure()
for nr, marker in enumerate(plc.regionMarker()):
    print(
        'Position marker number {}:'.format(nr + 1),
        marker.x(), marker.y(), marker.z()
    )
    ax.scatter(marker.x(), marker.y(), s=(4 - nr) * 20, color='k')
ax.set_title('marker positions')
fig.savefig('polygons.pdf')
###############################################################################
rect1 = mt.createRectangle(
    start=[0.0, 0.0],
    end=[2.0, -1.0],
    isClosed=True,
    marker=1,
)

# move the rectangle by changing the center position
rect2 = mt.createRectangle(
    start=[4.0, 0.0],
    end=[6.0, -1.0],
    isClosed=True,
    marker=2,
)

# scaled version
rect3 = mt.createRectangle(
    pos=[2.0, -2.75],
    size=[4.0, 2.5],
    isClosed=True,
    marker=3,
)

plc = mt.mergePLC([rect1, rect2, rect3])

ax, cb = pg.show(plc, markers=True)

fig = ax.get_figure()
for nr, marker in enumerate(plc.regionMarker()):
    print(
        'Position marker number {}:'.format(nr + 1),
        marker.x(), marker.y(), marker.z()
    )
    ax.scatter(marker.x(), marker.y(), s=(4 - nr) * 20, color='k')
ax.set_title('marker positions')
fig.savefig('rect.pdf')

###############################################################################
rect1 = mt.createRectangle(
    start=[0.0, 0.0],
    end=[2.0, -1.0],
    isClosed=True,
    marker=1,
)

# move the rectangle by changing the center position
rect2 = mt.createRectangle(
    start=[0.0, 0.0],
    end=[1.0, -0.5],
    isClosed=True,
    marker=2,
)

plc = mt.mergePLC([rect1, rect2])

ax, cb = pg.show(plc, markers=True)

fig = ax.get_figure()
for nr, marker in enumerate(plc.regionMarker()):
    print(
        'Position marker number {}:'.format(nr + 1),
        marker.x(), marker.y(), marker.z()
    )
    ax.scatter(marker.x(), marker.y(), s=(4 - nr) * 20, color='k')
ax.set_title('marker positions')
fig.savefig('rect_bad.pdf')

###############################################################################
rect1 = mt.createRectangle(
    start=[0.0, 0.0],
    end=[2.0, -1.0],
    isClosed=True,
    marker=1,
    markerPosition=[1.75, -0.25],
)

# move the rectangle by changing the center position
rect2 = mt.createRectangle(
    start=[0.0, 0.0],
    end=[1.0, -0.5],
    isClosed=True,
    marker=2,
)

plc = mt.mergePLC([rect1, rect2])

ax, cb = pg.show(plc, markers=True)

fig = ax.get_figure()
for nr, marker in enumerate(plc.regionMarker()):
    print(
        'Position marker number {}:'.format(nr + 1),
        marker.x(), marker.y(), marker.z()
    )
    ax.scatter(marker.x(), marker.y(), s=(4 - nr) * 20, color='k')
ax.set_title('marker positions')
fig.savefig('rect_fixed.pdf')


###############################################################################
# .. note::
#
#   This example was kindly contributed by Maximilian Weigand (University of
#   Bonn). If you also want to contribute an interesting example, check out
#   our `contribution guidelines <https://www.pygimli.org/contrib.html>`_.
