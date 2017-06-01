#!/usr/bin/env python

"""
    Mesh extrusion from 1D over 2D to a 3D hex mesh
"""

import math
import numpy as np
import pygimli as pg
from pygimli.meshtools import polytools as plc
# from pygimli.viewer import *
# import pygimli.meshtools as mt


c1 = plc.createCircle(pos=(-5.0, 0.0), radius=0.1, segments=5,
                      start=math.pi, end=2*math.pi, isClosed=False)
c2 = plc.createCircle(pos=(5.0, 0.0), radius=0.1, segments=5,
                      start=math.pi, end=2*math.pi, isClosed=False)

left = plc.createLine(start=(-20, 0.0), end=(-5.1, 0.0), segments=10)
left.node(8).setMarker(1)
mid = plc.createLine(start=(-4.9, 0.0), end=(4.9, 0.0), segments=20)
right = plc.createLine(start=(5.1, 0.0), end=(20, 0.0), segments=10)
left.node(2).setMarker(1)

border = plc.mergePLC([left, c1, mid, c2, right])

depth = 20
nz = 15
newNodes = []
y = pg.increasingRange(0.2, depth, nz)
surface = pg.createMesh2D(border, y, 0, 0, False)

#for n in surface.nodes():
    #yNodes = pg.increasingRange(yDefault[1], depth+n.y(), nz)
    #for y in yNodes[0:]:
        #newNodes.append([n.x(), n.y() -y])

#surface = pg.createGrid(x=yDefault, y=pg.sort(pg.x(surface.positions())))

#for i, n in enumerate(surface.nodes()):
    #n.setPos(newNodes[i])

#surface.smooth(1, 1, 1, 10)

ax, _ = pg.show(surface)
#showBoundaryNorm(surface, axes=ax)

ax.set_ylim([-21, 1])
ax.set_xlim([-21, 21])

zinc = pg.increasingRange(1, 10, 5)
zmid = np.linspace(1, 10, 20)
z = pg.cat(-(np.asarray(zinc)[::-1]), zmid)
z = pg.cat(z, 11+zinc)

mesh3D = pg.createMesh3D(surface, z)
mesh3D.exportVTK("mesh3d")

pg.wait()
