import numpy as np
import pygimli as pg
from .modelling import MeshModelling


class MultiFrameModelling(MeshModelling):
    """Full frame (multiple fop parallel) forward modelling."""
    def __init__(self, scalef=1.0):
        """Init class and jacobian matrix."""
        super().__init__()
        self.jac = pg.matrix.BlockMatrix()
        self.scalef = scalef

    def setData(self, alldata, fop):
        """Distribute the data containers amongst the fops."""
        self.fops = []
        for i, data in enumerate(alldata):
            fopi = fop()
            fopi.setData(data)
            self.fops.append(fopi)

    def setMeshPost(self, mesh):
        for i, fop in enumerate(self.fops):
            fop.setMesh(mesh, ignoreRegionManager=True)

        self.prepareJacobian()

    def setDefaultBackground(self):
        """Set the default background behaviour."""
        regionIds = self.regionManager().regionIdxs()
        pg.info("Found {} regions.".format(len(regionIds)))
        if len(regionIds) > 1:
            bk = pg.sort(regionIds)[0]
            pg.info("Region with smallest marker ({0}) "
                    "set to background".format(bk))
            self.setRegionProperties(bk, background=True)

    @property
    def parameterCount(self):
        return self.regionManager().parameterCount() * len(self.fops)

    def prepareJacobian(self):
        """Build up Jacobian block matrix (once the sizes are known)."""
        self.jac.clear()
        self.nm = self.regionManager().parameterCount()
        self.nf = len(self.fops)
        print(self.nm, "model cells")
        nd = 0
        for i, fop in enumerate(self.fops):
            self.jac.addMatrix(fop.jacobian(), nd, i*self.nm)
            nd += fop.data.size()

        self.jac.recalcMatrixSize()
        self.setJacobian(self.jac)

    def createConstraints(self):
        """Create constraint matrix (special type for this)."""
        if isinstance(super().createConstraints(), pg.SparseMapMatrix):
            self.C1 = pg.SparseMapMatrix(self.constraintsRef())
            # make a copy because it will be overwritten
        else:
            self.C1 = self.constraints()

        self.C = pg.matrix.FrameConstraintMatrix(self.C1,
                                                 len(self.fops),
                                                 self.scalef)
        self.setConstraints(self.C)
        # cw = self.regionManager().constraintWeights()
        # self.regionManager().setConstraintsWeights(np.tile(cw, self.nf))
        ## switch off automagic inside core.inversion which checks for local modeltransform of the regionManager
        self.regionManager().setLocalTransFlag(False)

    def response(self, model):
        mod = np.reshape(model, [len(self.fops), -1])
        return np.concatenate([fop.response(mo) for fop, mo in
                               zip(self.fops, mod)])

    def createJacobian(self, model):
        mod = np.reshape(model, [len(self.fops), -1])
        for i, fop in enumerate(self.fops):
            fop.createJacobian(mod[i])

    def createDefaultStartModel(self):  # , dataVals):
        return pg.Vector(self.nm*self.nf, 10.0)  # look up in fop
    #     return np.concatenate([fop.createDefaultStartModel() for fop in self.fops])

    def createStartModel(self, dataVals):
        # return np.concatenate([fop.createStartModel() for fop in self.fops])
        self.nm = self.regionManager().parameterCount()
        return pg.Vector(self.nm*len(self.fops), np.median(dataVals))
