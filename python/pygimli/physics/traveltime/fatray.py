import numpy as np
import pygimli as pg


class FatrayDijkstraModelling(pg.TravelTimeDijkstraModelling):
    """Shortest-path (Dijkstra) based travel time with fat ray jacobian."""
    def __init__(self, frequency=1000., verbose=False):
        super().__init__(verbose)
        self.frequency = frequency

    def setMesh(self, mesh, secondaryNodes=3):
        """Set mesh and create secondary Nodes."""
        super().setMesh(mesh)
#        self.mesh0 = mesh
#        self.mesh1 = self.mesh().createSecondaryNodes(secondaryNodes)
#        self.setMesh(self.mesh1, ignoreRegionManager=True)
        self.IM = self.mesh().interpolationMatrix(self.mesh().cellCenters())

    def createJacobian(self, slowness):
        """Generate Jacobian matrix using fat-ray after Jordi et al. (2016)."""
        self.J = pg.Matrix(self.data().size(), self.mesh().cellCount())
        self.sensorNodes = [self.mesh().findNearestNode(pos)
                            for pos in self.data().sensorPositions()]
        Di = self.dijkstra()
        numN = self.mesh().nodeCount()
        data = self.data()
        numS = data.sensorCount()
        Tmat = pg.RMatrix(numS, numN)
        Dmat = pg.RMatrix(numS, numS)
        for i, node in enumerate(self.sensorNodes):
            Di.setStartNode(node)
            Tmat[i] = Di.distances()(0, numN)
            Dmat[i] = Tmat[i][self.sensorNodes]

        for i in range(data.size()):
            iS = int(data("s")[i])
            iG = int(data("g")[i])
            tsr = Dmat[iS][iG]  # shot-receiver travel time
            dt = self.IM * (Tmat[iS] + Tmat[iG]) - tsr
            weight = np.maximum(1 - 2 * self.frequency * dt, 0.0)  # 1 on ray
            wa = weight  # * np.sqrt(self.mesh().cellSizes())
            if np.sum(wa) > 0:  # not if all values are zero
                wa /= np.sum(wa)

            self.J[i] = wa * tsr / slowness

        self.setJacobian(self.J)
