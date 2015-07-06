#!/usr/bin/env python
"""
    Fast Marching test
"""
import pygimli as pg
import numpy as np
import matplotlib.pyplot as plt
import time
from pygimli.mplviewer import drawMesh, drawField, drawStreamLines
#import heapq
from math import asin, tan

from fastMarchingTest import fastMarch

"""
Solve the particular Hamilton-Jacobi (HJ) equation, known as the Eikonal equation
.. math::
    |\grad u(x)| & = f(x) \\
    ||\grad t||_2 &= s
    \grad t^2 &= s^2~\cite{SunFomel2000}

where :math:`t` denote traveltime for a spatial distributed slowness :math:`s`

In the special case when f(x) = 1, the solution gives the signed distance from the boundary
"""

# def fastMarch(mesh, downwind, times, upTags, downTags):

# def findSlowness(edge):
# if edge.leftCell() is None: # outer boundary
#slowness = edge.rightCell().attribute()
# elif edge.rightCell() is None: # does not
#slowness = edge.leftCell().attribute()
# else: # take fastest (smallest slowness)
#slowness = min(edge.leftCell().attribute(), edge.rightCell().attribute())
# return slowness
# def findSlowness(...)

#upCandidate = []

# for node in downwind:
#neighNodes = pg.commonNodes(node.cellSet())

#upNodes = []
# for n in neighNodes:
# if upTags[ n.id() ]:
# upNodes.append(n)

# if len(upNodes) == 1:
# this is the dijkstra case
#edge = pg.findBoundary(upNodes[0], node)
#tt = times[ upNodes[0].id() ] + findSlowness(edge) * edge.shape().domainSize()

#heapq.heappush(upCandidate, (tt, node))
# else:
#cells = node.cellSet()
# for c in cells:
# for i in range(c.nodeCount()):
#edge = pg.findBoundary(c.node(i), c.node((i + 1)%3))

#a = edge.node(0)
#b = edge.node(1)
#ta = times[ a.id() ]
#tb = times[ b.id() ]

# if upTags[ a.id() ] and upTags[ b.id() ]:
#line = pg.Line(a.pos(), b.pos())
#t = min(1., max(0., line.nearest(node.pos())))

#ea = pg.findBoundary(a, node)
#eb = pg.findBoundary(b, node)

# if t == 0:
#slowness = findSlowness(ea)
# elif t == 1:
#slowness = findSlowness(eb)
# else:
#slowness = c.attribute()

#ttimeA = (ta               + slowness * a.pos().distance(node.pos()))
#ttimeQ = (ta + t*(tb-ta)) + slowness * line(t).distance(node.pos())
#ttimeB = (tb               + slowness * b.pos().distance(node.pos()))

#heapq.heappush(upCandidate, (min(ttimeA,ttimeQ,ttimeB), node))

#candidate = heapq.heappop(upCandidate)
#newUpNode = candidate[1]
#times[ newUpNode.id() ] = candidate[0]
#upTags[ newUpNode.id() ] = 1
# downwind.remove(newUpNode)

#newDownNodes = pg.commonNodes(newUpNode.cellSet())
# for nn in newDownNodes:
# if not upTags[ nn.id() ] and not downTags[ nn.id() ]:
# downwind.add(nn)
#downTags[ nn.id() ] = 1

# def fastMarch(...)

if __name__ == '__main__':
    xmin, xmax, zlay = -20., 150., 20.
    # create PLC (piece-wise linear complex) of two layers
    PLC = pg.Mesh(2)
    nodes = []
    #   0-----------1
    #   |           |
    #   5-----------2
    #   |           |
    #   4-----------3
    nodes.append(PLC.createNode(xmin, 0., 0.))
    nodes.append(PLC.createNode(xmax, 0., 0.))
    nodes.append(PLC.createNode(xmax, -zlay, 0.))
    nodes.append(PLC.createNode(xmax, -zlay * 2, 0.))
    nodes.append(PLC.createNode(xmin, -zlay * 2, 0.))
    nodes.append(PLC.createNode(xmin, -zlay, 0.))
    # connect the nodes
    for i in range(5):
        PLC.createEdge(nodes[i], nodes[i + 1])

    PLC.createEdge(nodes[5], nodes[0])
    PLC.createEdge(nodes[5], nodes[2])
    PLC.addRegionMarker(pg.RVector3(0., -zlay + .1), 0, 10.)
    PLC.addRegionMarker(pg.RVector3(0., -zlay - .1), 1, 10.)
    # insert region markers into the two layers and make mesh
    tri = pg.TriangleWrapper(PLC)
    tri.setSwitches('-pzeAfaq34.5')
    mesh = pg.Mesh(2)
    tri.generate(mesh)
    mesh.createNeighbourInfos()
    print(mesh)

    # make velocity model
    v = [1000., 3000.]
    slomap = pg.stdMapF_F()  # map for mapping real slowness values
    for i, vi in enumerate(v):
        slomap.insert(i, 1. / vi)

    mesh.mapCellAttributes(slomap)  # map values to attributes using map

    # initialize source position and trvel time vector
    source = pg.RVector3(0., 0.)
    times = pg.RVector(mesh.nodeCount(), 0.)

    # initialize sets and tags
    upwind, downwind = set(), set()
    upTags, downTags = np.zeros(mesh.nodeCount()), np.zeros(mesh.nodeCount())

    # define initial condition
    cell = mesh.findCell(source)

    for i, n in enumerate(cell.nodes()):
        times[n.id()] = cell.attribute() * n.pos().distance(source)
        upTags[n.id()] = 1
    for i, n in enumerate(cell.nodes()):
        tmpNodes = pg.commonNodes(n.cellSet())
        for nn in tmpNodes:
            if not upTags[nn.id()] and not downTags[nn.id()]:
                downwind.add(nn)
                downTags[nn.id()] = 1

    # start fast marching
    tic = time.time()
    while len(downwind) > 0:
        fastMarch(mesh, downwind, times, upTags, downTags)

    print(time.time() - tic, "s")

    # compare with analytical solution along the x axis
    x = np.arange(-20., 150., 0.5)
    t = pg.interpolate(mesh, times, pg.asvector(x), x * 0., x * 0.)
    tdirect = np.abs(x) / v[0]  # direct wave
    alfa = asin(v[0] / v[1])  # critically refracted wave angle
    xreflec = tan(alfa) * zlay * 2.  # first critically refracted
    trefrac = (np.abs(x) - xreflec) / v[1] + xreflec * v[1] / v[0]**2
    tana = np.where(trefrac < tdirect, trefrac, tdirect)  # minimum of both
    print("min(dt)=", min(t-tana)*1e3, "ms max(dt)=", max(t-tana)*1e3, "ms")

    # %% plot traveltime field, a few lines
    fig, ax = plt.subplots(nrows=2, sharex=True)
    drawMesh(ax[0], mesh)
    drawField(ax[0], mesh, times, True)
    drawStreamLines(ax[0], mesh, times, color='white')
    # plot calculated and measured travel times
    ax[1].plot(x, t, 'b.-', x, tdirect, 'r-', x, trefrac, 'g-')
    plt.show()
