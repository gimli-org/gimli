#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Tue Jul 07 09:44:36 2015

@author: Marcus
"""

import heapq
import pygimli as pg
import numpy as np
from pygimli.physics.traveltime.ratools import createGradientModel2D


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
    upCandidate = []
#    print('.', end='')

    for node in downwind:
        neighNodes = pg.commonNodes(node.cellSet())

        upNodes = []
        for n in neighNodes:
            if upT[n.id()]:
                upNodes.append(n)

        if len(upNodes) == 1:  # the Dijkstra case
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

    candidate = heapq.heappop(upCandidate)
    newUpNode = candidate[1]  # original
    times[newUpNode.id()] = candidate[0]
    upT[newUpNode.id()] = 1
    downwind.remove(newUpNode)
    newDownNodes = pg.commonNodes(newUpNode.cellSet())
#    newUpNodeId = candidate[1]  # original
#    times[newUpNodeId] = candidate[0]
#    upT[newUpNodeId] = 1
#    downwind.remove(newUpNodeId)
#    newDownNodes = pg.commonNodes(mesh.node(newUpNodeId).cellSet())
    for nn in newDownNodes:
        if not upT[nn.id()] and not downT[nn.id()]:
            downwind.add(nn)
            downT[nn.id()] = 1


class TravelTimeFMM(pg.ModellingBase):
    """Modelling class using the Fast Marching Method (FMM).

    It can be used alternatively to Dijkstra modelling.
    However, currently it is quite slow.
    A implementation in C++ might speed up.
    """
    def __init__(self, mesh=None, data=None, frequency=200, verbose=False):
        """
        Init function.

        Parameters:
        -----------
        mesh : pygimli.Mesh
            2D mesh to be used in the forward calculations.
        data : pygimli.DataContainer
            The datacontainer with sensor positions etc.
        verbose : boolean
            More printouts or not...
        """

        pg.ModellingBase.__init__(self, verbose)
        super().__init__(verbose=verbose)
        self.debug = False
        if mesh is not None:
            self.setMesh(mesh)  # besser use createRefinedForwardMesh
        if data is None:
            self.setData(pg.DataContainer())
        else:
            self.setData(data)
        self.frequency = frequency

    def setMesh(self):
        """Set mesh"""
        super().setMesh(mesh)
        self.timeMatrix = np.zeros((self.data().sensorCount(),
                                    self.mesh().cellCount()))

    def setData(self, data):
        """Set data and prepare stuff"""
        super().setData(data)
        self.nModel = mesh.cellCount()
        self.midPoints = mesh.cellCenters()
        self.nNodes = self.mesh().nodeCount()
        self.cellSizes = self.mesh_.cellSizes()
        nSensors = self.data.sensorCount()
        nModels = self.mesh().cellCount()
        self.dataMatrix = np.zeros((nSensors, nSensors))
        self.timeMatrix = np.zeros((nSensors, nModels))

    def computeTravelTimes(self, slowness, calcOthers=False):
        """Compute the travel times and fill data and time matrix
        for later use of response and Jacobian, respectively.
        For response only active sources are needed, for Jacobian all."""
        # mesh = self.mesh()  # better but for now input mesh
        mesh = self.mesh_
        param_markers = np.unique(mesh.cellMarkers())
        param_count = len(param_markers)
        if len(slowness) == mesh.cellCount():
            mesh.setCellAttributes(slowness)
            # self.mapModel(slowness)
        elif len(slowness) == param_count:
            # map the regions in the mesh to slowness
            slow_map = pg.stdMapF_F()
            min_reg_num = min(param_markers)
            for i, si in enumerate(slowness):
                slow_map.insert(float(i+min_reg_num), si)

            mesh.mapCellAttributes(slow_map)
        else:
            raise ValueError("Wrong no of parameters. Mesh size: {}, no "
                             "of regions: {}, and number of slowness values:"
                             "{}".format(self.mesh().cellCount(), param_count,
                                         len(slowness)))

        times = pg.RVector(self.nNodes, 0.)
        upTags = np.zeros(self.nNodes)
        downTags = np.zeros(mesh.nodeCount())
        sourceIndices = np.unique(self.data_("s"))
        if calcOthers:
            ns = len(sourceIndices)
            geophoneIndices = np.setxor1d(np.arange(self.data_.sensorCount()),
                                          sourceIndices)
            sourceIndices = geophoneIndices
#            geophoneIndices = np.unique(self.data_("g"))
            print("{:d}-{:d}={:d}".format(
                self.data_.sensorCount(), ns, len(sourceIndices)))
        if self.debug:  # resize not working
            self.solution().resize(self.mesh().nodeCount(), self.nSensors)
            print(self.solution().rows(), self.solution().cols())
        for iSource in np.array(sourceIndices, dtype=int):
            if self.debug:
                print(iSource)
            # initial condition (reset vectors)
            times *= 0.0
            upTags *= 0
            downTags *= 0
            downwind = set()
            source = self.data_.sensorPosition(int(iSource))
            cell = self.mesh_.findCell(source)
            # fill in nodes around source using local smoothness
            for i, n in enumerate(cell.nodes()):
                times[n.id()] = cell.attribute() * n.pos().distance(source)
                upTags[n.id()] = 1
            for i, n in enumerate(cell.nodes()):
                tmpNodes = pg.commonNodes(n.cellSet())
                for nn in tmpNodes:
                    if not upTags[nn.id()] and not downTags[nn.id()]:
                        downwind.add(nn)
                        downTags[nn.id()] = 1

            while len(downwind) > 0:  # start fast marching
                fastMarch(self.mesh_, downwind, times, upTags, downTags)

            self.dataMatrix[iSource] = pg.interpolate(
                mesh, times, self.data_.sensorPositions())
            self.timeMatrix[iSource] = pg.interpolate(
                mesh, times, self.midPoints)
            if self.debug:
                print(self.solution().rows(), self.solution().cols())
                print(len(times), self.mesh())
                self.solution()[int(iSource)] = times
                self.solution().setCol(int(iSource), times)

    def response(self, slowness):
        """
        Response function. Returns the result of the forward calculation.
        Uses the shot- and sensor positions specified in the data container.
        """
        self.computeTravelTimes(slowness)
        # assembling the data from the data matrix
        data = self.data_
        n_data = data.size()
        t_fmm = pg.RVector(n_data)
        for i in range(n_data):
            t_fmm[i] = self.dataMatrix[int(data("s")[i])][int(data("g")[i])]

        return t_fmm

    def createJacobian(self, slowness):
        """
        Jacobian matrix using a fat-ray approach (Jordi et al. 2016).
        """
        self.jacobian().resize(self.data.size(), self.mesh().cellCount())
        # first compute reciprocal travel times for geophone sources
        self.computeTravelTimes(slowness, calcOthers=True)
        data = self.data_
        n_data = data.size()
#        slo = mesh.cellAttributes()  # in case of regions
        for i in range(n_data):
            iS, iG = int(data("s")[i]), int(data("g")[i])
            tsr = self.dataMatrix[iS][iG]
            dt = self.timeMatrix[iS] + self.timeMatrix[iG] - tsr
            weight = np.maximum(1 - 2 * self.frequency * dt, 0.0)
#            print(pg.sum(pg.sign(weight)))
            wa = weight * self.cellSizes
            self.jacobian()[i] = wa / np.sum(wa) * tsr / slowness

if __name__ == '__main__':
    """
    Set up FMM modelling operator and run a synthetic model
    """
    data = pg.DataContainer('example_topo.sgt', 's g')
    print(data)
    mesh = pg.meshtools.createParaMesh(data, boundary=0, paraBoundary=5,
                                       paraDepth=20,
                                       quality=34.5, paraMaxCellSize=5)
    mesh.createNeighbourInfos()
    print(mesh)
    slo = createGradientModel2D(data, mesh, vTop=500, vBot=2500)
    fwd = TravelTimeFMM(mesh, data, frequency=500)  #
    fwd.createRefinedForwardMesh(False)
    resp = fwd.response(slo)
    data.set('t', resp)
    print("ready with response, starting jacobian")
    fwd.createJacobian(slo)
    raise SystemExit
    # %%
    pg.plt.imshow(fwd.dataMatrix, interpolation='nearest')
    fwd.createJacobian(slo)
    one = pg.RVector(data.size(), 1.0)
    coverage = fwd.jacobian().transMult(one)
