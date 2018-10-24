import numpy as np
import pygimli as pg


class FatrayDijkstraModellingInterpolate(pg.TravelTimeDijkstraModelling):
    """Shortest-path (Dijkstra) based travel time with fat ray jacobian."""
    def __init__(self, frequency=100., verbose=False):
        super().__init__(verbose)
        self.frequency = frequency
        self.iMat = pg.SparseMapMatrix()

    def createJacobian(self, slowness):
        """Generate Jacobian matrix using fat-ray after Jordi et al. (2016)."""
        self.J = pg.Matrix(self.data().size(), self.mesh().cellCount())
        self.sensorNodes = [self.mesh().findNearestNode(pos)
                            for pos in self.data().sensorPositions()]
        if (self.iMat.cols() != self.mesh().nodeCount() or
            self.iMat.rows() != self.mesh().cellCount()):
            self.iMat = self.mesh().interpolationMatrix(
                    self.mesh().cellCenters())
        Di = self.dijkstra()
        slowPerCell = self.createMappedModel(slowness, 1e16)
        Di.setGraph(self.createGraph(slowPerCell))
        numN = self.mesh().nodeCount()
        data = self.data()
        numS = data.sensorCount()
        Tmat = pg.RMatrix(numS, numN)
        Dmat = pg.RMatrix(numS, numS)
        for i, node in enumerate(self.sensorNodes):
            Di.setStartNode(node)
            Tmat[i] = Di.distances()  # (0, numN)
            Dmat[i] = Tmat[i][self.sensorNodes]

        for i in range(data.size()):
            iS = int(data("s")[i])
            iG = int(data("g")[i])
            tsr = Dmat[iS][iG]  # shot-receiver travel time
            dt = self.iMat * (Tmat[iS] + Tmat[iG]) - tsr
            weight = np.maximum(1 - 2 * self.frequency * dt, 0.0)  # 1 on ray
            wa = weight  # * np.sqrt(self.mesh().cellSizes())
            if np.sum(wa) > 0:  # not if all values are zero
                wa /= np.sum(wa)

            self.J[i] = wa * tsr / slowness

        self.setJacobian(self.J)


class FatrayDijkstraModellingMidpoint(pg.TravelTimeDijkstraModelling):
    """Shortest-path (Dijkstra) based travel time with fat ray jacobian."""
    def __init__(self, frequency=100., verbose=False):
        super().__init__(verbose)
        self.frequency = frequency

    def setMesh(self, mesh, **kwargs):  # secondaryNodes=3):
        """Set mesh and create additional secondary Nodes in cell centers."""
        super().setMesh(mesh, **kwargs)  # ignoreRegionManager=True)
        print(self.mesh(), self.mesh().secondaryNodeCount())
        self.mids = pg.IVector()
        self.nnodes = self.mesh().nodeCount()
        for c in self.mesh().cells():
            n = self.mesh().createSecondaryNode(c.center())
            c.addSecondaryNode(n)
            self.mids.push_back(n.id())

        print(self.mesh())

    def createJacobian(self, slowness):
        """Generate Jacobian matrix using fat-ray after Jordi et al. (2016)."""
        self.J = pg.Matrix(self.data().size(), self.mesh().cellCount())
        self.sensorNodes = [self.mesh().findNearestNode(pos)
                            for pos in self.data().sensorPositions()]
        Di = self.dijkstra()
        slowPerCell = self.createMappedModel(slowness, 1e16)
        Di.setGraph(self.createGraph(slowPerCell))
        numN = self.mesh().nodeCount()
        data = self.data()
        numS = data.sensorCount()
        Tmat = pg.RMatrix(numS, numN)
        Dmat = pg.RMatrix(numS, numS)
        print(self.mesh())
        print(self.nnodes, max(self.mids))
        for i, node in enumerate(self.sensorNodes):
            Di.setStartNode(node)
            dist0 = Di.distances()
            dist = Di.distances(withSecNodes=True)
            print("dist len ", len(dist0), len(dist))
            Tmat[i] = dist[self.mids]
#            Tmat[i] = (self.nnodes, len(dist))
            Dmat[i] = Tmat[i][self.sensorNodes]

        for i in range(data.size()):
            iS = int(data("s")[i])
            iG = int(data("g")[i])
            tsr = Dmat[iS][iG]  # shot-receiver travel time
            dt = Tmat[iS] + Tmat[iG] - tsr
            weight = np.maximum(1 - 2 * self.frequency * dt, 0.0)  # 1 on ray
            wa = weight  # * np.sqrt(self.mesh().cellSizes())
            if np.sum(wa) > 0:  # not if all values are zero
                wa /= np.sum(wa)

            self.J[i] = wa * tsr / slowness

        self.setJacobian(self.J)


FatrayDijkstraModelling = FatrayDijkstraModellingInterpolate
# FatrayDijkstraModelling = FatrayDijkstraModellingMidpoint
