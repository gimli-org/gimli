#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
Solve the particular Hamilton-Jacobi (HJ) equation (Eikonal equation).

Solve the particular Hamilton-Jacobi (HJ) equation, known as the Eikonal
equation :cite:`SunFomel2000`

.. math::

    |\grad u(x)| & = f(x) \\\\
    |\grad t| & = s \\\\
    (\grad t)^2 & = s^2

.. error::
    check equation and reference!!!

where :math:`t` denote traveltime for a spatial distributed slowness :math:`s`

In the special case when :math:`f(x) = 1`, the solution gives the signed
distance from the boundary
"""

import time
import heapq

import numpy as np
import matplotlib.pyplot as plt

import pygimli as pg
from pygimli.mplviewer import drawMesh, drawField


def findSlowness(edge):
    """WRITEME."""
    if edge.leftCell() is None:
        slowness = edge.rightCell().attribute()
    elif edge.rightCell() is None:
        slowness = edge.leftCell().attribute()
    else:
        slowness = min(edge.leftCell().attribute(),
                       edge.rightCell().attribute())
    return slowness
# def findSlowness(...)


def fastMarch(mesh, downwind, times, upT, downT):
    """WRITEME."""
    print(mesh)
    upCandidate = []

    for node in downwind:
        neighNodes = pg.commonNodes(node.cellSet())

        upNodes = []
        for n in neighNodes:
            if upT[n.id()]:
                upNodes.append(n)

        if len(upNodes) == 1:
            # this is the dijkstra case
            edge = pg.findBoundary(upNodes[0], node)
            tt = times[upNodes[0].id()] + \
                findSlowness(edge) * edge.shape().domainSize()

            heapq.heappush(upCandidate, (tt, node))
        else:
            cells = node.cellSet()
            for c in cells:
                for i in range(c.nodeCount()):
                    edge = pg.findBoundary(c.node(i), c.node((i + 1) % 3))

                    a = edge.node(0)
                    b = edge.node(1)
                    ta = times[a.id()]
                    tb = times[b.id()]

                    if upT[a.id()] and upT[b.id()]:
                        line = pg.Line(a.pos(), b.pos())
                        t = min(1., max(0., line.nearest(node.pos())))

                        ea = pg.findBoundary(a, node)
                        eb = pg.findBoundary(b, node)

                        if t == 0:
                            slowness = findSlowness(ea)
                        elif t == 1:
                            slowness = findSlowness(eb)
                        else:
                            slowness = c.attribute()

                        ttimeA = (ta + slowness * a.pos().distance(node.pos()))
                        ttimeQ = (ta + t * (tb - ta)) + \
                            slowness * line(t).distance(node.pos())
                        ttimeB = (tb + slowness * b.pos().distance(node.pos()))

                        heapq.heappush(upCandidate,
                                       (min(ttimeA, ttimeQ, ttimeB), node))

    # for c in upCandidate:
        # print c[1].id(), c[0]
    candidate = heapq.heappop(upCandidate)
    # print candidate
    newUpNode = candidate[1]
    times[newUpNode.id()] = candidate[0]
    upT[newUpNode.id()] = 1
    # print newUpNode
    downwind.remove(newUpNode)

    newDownNodes = pg.commonNodes(newUpNode.cellSet())
    for nn in newDownNodes:
        if not upT[nn.id()] and not downT[nn.id()]:
            downwind.add(nn)
            downT[nn.id()] = 1


def test():
    """WRITEME."""
    mesh = pg.Mesh('mesh/test2d')
    mesh.createNeighbourInfos()

    print(mesh)

    source = pg.RVector3(-80, 0.)
    times = pg.RVector(mesh.nodeCount(), 0.)

    for c in mesh.cells():
        if c.marker() == 1:
            c.setAttribute(1.)
        elif c.marker() == 2:
            c.setAttribute(0.5)
        # c.setAttribute(abs(1./c.center()[1]))

    fig = plt.figure()
    ax = fig.subplot(1, 1, 1)

    anaTimes = pg.RVector(mesh.nodeCount(), 0.0)

    for n in mesh.nodes():
        anaTimes[n.id()] = source.distance(n.pos())

    # d = pg.DataContainer()
    # dijk = pg.TravelTimeDijkstraModelling(mesh, d)

    # upwind = set()
    downwind = set()
    upTags = np.zeros(mesh.nodeCount())
    downTags = np.zeros(mesh.nodeCount())

    # define initial condition
    cell = mesh.findCell(source)

    for n in cell.nodes():
        times[n.id()] = cell.attribute() * n.pos().distance(source)
        upTags[n.id()] = 1

    for n in cell.nodes():
        tmpNodes = pg.commonNodes(n.cellSet())
        for nn in tmpNodes:
            if not upTags[nn.id()] and not downTags[nn.id()]:
                downwind.add(nn)
                downTags[nn.id()] = 1

    # start fast marching
    tic = time.time()
    while len(downwind) > 0:
        # model = pg.RVector(mesh.cellCount(), 0.0)

        # for c in upwind: model.setVal(2, c.id())
        # for c in downwind: model.setVal(1, c.id())
        # drawMesh(a, mesh)

        # a.plot(source[0], source[1], 'x')

        # for c in mesh.cells():
        # a.text(c.center()[0],c.center()[1], str(c.id()))

        # for n in mesh.nodes():
        # if upTags[n.id()]:
        # a.plot(n.pos()[0], n.pos()[1], 'o', color='black')

        # mindist = 9e99
        # for n in downwind:
        # a.plot(n.pos()[0], n.pos()[1], 'o', color='white')

        # if n.pos().distance(source) < mindist:
        # mindist = n.pos().distance(source)
        # nextPredict = n

        # print "next predicted ", n.id(), mindist
        # a.plot(nextPredict.pos()[0], nextPredict.pos()[1],
        # 'o', color='green')

        # drawField(a, mesh, times)
        # drawField(a, mesh, anaTimes, colors = 'white')
        # a.figure.canvas.draw()

        # a.figure.show()
        # raw_input('Press Enter...')
        # a.clear()

        fastMarch(mesh, downwind, times, upTags, downTags)

    print(time.time() - tic, "s")

    drawMesh(ax, mesh)
    drawField(ax, mesh, times, filled=True)
    # drawStreamCircular(a, mesh, times, source, 30.,
    #                nLines = 50, step = 0.1, showStartPos = True)

    # ax1.streamplot(X, Y, U, V, density=[0.5, 1])

    # drawStreamLinear(a, mesh, times,
    # pg.RVector3(-100., -10.0),
    # pg.RVector3(100., -10.0),
    # nLines = 50, step = 0.01,
    # showStartPos = True)


if __name__ == '__main__':
    test()
    plt.show()
