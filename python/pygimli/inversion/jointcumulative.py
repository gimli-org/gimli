# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 16:58:13 2015

@author: Guenther.T
"""

import pygimli as pg
import pybert as pb


class RMultRMatrix(pg.MatrixBase):  # shift to C++ or python/matrices
    """ matrix A to be multiplied by a right hand side vector r """
    def __init__(self, A, r):
        super().__init__()
        self.A = A
        self.r = r

    def mult(self, x):
        """ return M*x = A*(r*x) """
        return self.A.mult(x*self.r)

    def transMult(self, x):
        """ return A*x """
        return self.A.transMult(x) * self.r

    def cols(self):
        """ number of columns """
        return self.A.cols()

    def rows(self):
        """ number of rows """
        return self.A.rows()


class jointCrossGradientsModelling(pg.ModellingBase):
    """ cumulative forward operator """
    def __init__(self, mesh, dataERT, dataSRT, verbose=True):
        super().__init__(verbose=verbose)
        #  set up two individual forward operators
        self.fERT = pb.DCSRMultiElectrodeModelling(mesh, dataERT)
        self.fSRT = pg.TravelTimeDijkstraModelling(mesh, dataSRT)
        for fi in [self.fERT, self.fSRT]:
            fi.region(1).setBackground(True)
            fi.createRefinedForwardMesh(True)
            fi.createConstraints()
        # save inversion mesh and retrieve number of (single) unknowns
        self.pd = self.fERT.regionManager().paraDomain()
        self.nModel = self.pd.cellCount()
        print("model: " + str(self.nModel) + "*2=" + str(self.nModel*2))
        # as the fop has no mesh it needs to know the number of unknowns
        self.regionManager().setParameterCount(self.nModel*2)
        self.nBounds = self.fERT.constraints().nRows()
        # set up Jacobian as block-diagonal with two individual Jacobians
        self.J = pg.RBlockMatrix()
        nERT = self.J.addMatrix(self.fERT.jacobian())
        self.J.addMatrixEntry(nERT, 0, 0)
        nSRT = self.J.addMatrix(self.fSRT.jacobian())
        self.J.addMatrixEntry(nSRT, dataERT.size(), self.nModel)
        self.setJacobian(self.J)
        # set up constraint matrix by two single C1 constraints and a matrix
        self.C = pg.RBlockMatrix()
        cERT = self.C.addMatrix(self.fERT.constraints())
        self.C.addMatrixEntry(cERT, 0, 0)
        cSRT = self.C.addMatrix(self.fSRT.constraints())
        self.C.addMatrixEntry(cSRT, self.nBounds, 0)
        tmpModel = pg.RVector(self.nModel*2)  # for now a zero model
        self.A = pg.IdentityMatrix(self.nModel*2)  # replace by CG matrix
        self.B = RMultRMatrix(self.A, tmpModel)
        aCG = self.C.addMatrix(self.B)
        self.C.addMatrixEntry(aCG, self.nBounds*2, 0)
        self.C.recalcMatrixSize()
        print("C: ", self.C.rows(), self.C.cols())
        self.setConstraints(self.C)  # not yet working as C is still SparseMap
        # set up (unit) constraint and model weight vectors
        self.wm = pg.RVector(self.nModel*2, 1.0)  # model weight
        self.wc = pg.RVector(self.nBounds*2+self.nModel*2, 1.0)  # constr. w.

    def response(self, model):
        """ paste responses of both methods (and update rhs vector) """
        # self.B.r = model
        return pg.cat(self.fERT.response(model(0, self.nModel)),
                      self.fSRT.response(model(self.nModel, self.nModel*2)))

    def createJacobian(self, model):
        """ fill the individual jacobian matrices"""
        self.fERT.createJacobian(model(0, self.nModel))
        self.fSRT.createJacobian(model(self.nModel, self.nModel*2))
        self.B.r = model  # update CG matrix
        self.J.recalcMatrixSize()  # just to check
        print("J: ", self.J.rows(), self.J.cols())

# %% load data abd add errors
dataERT = pb.DataContainerERT('example.data')
dataSRT = pg.DataContainer('example.sgt', 's g')
print("data: " + str(dataERT.size()) + "+" + str(dataSRT.size()) + "=" +
      str(dataERT.size() + dataSRT.size()))
if 'err' not in dataERT.dataMap().keys():
    dataERT.set('err', pg.RVector(dataERT.size(), 0.02))  # relative
if 'err' not in dataSRT.dataMap().keys():
    dataSRT.set('err', 0.001 / dataSRT('t') + 0.01)  # downgrade near offsets
# %% create mesh, setupforward operator and inversion
mesh = pg.meshtools.createParaMesh(dataSRT, paraDX=0.5)
f = jointCrossGradientsModelling(mesh, dataERT, dataSRT)
INV = pg.RInversion(pg.cat(dataERT('rhoa'), dataSRT('t')), f, True)
INV.setRelativeError(pg.cat(dataERT('err'), dataSRT('err')))
startModel = pg.cat(pg.RVector(f.nModel, pg.median(dataERT('rhoa'))),
                    pg.RVector(f.nModel, 0.001))
INV.setModel(startModel)
tRes, tSlo = pg.RTransLogLU(1., 1000.), pg.RTransLogLU(0.0001, 0.002)
transModel = pg.RTransCumulative()
transModel.add(tRes, f.nModel)
transModel.add(tSlo, f.nModel)
INV.setTransModel(transModel)
INV.setMWeight(f.wm)
INV.setCWeight(f.wc)
# %% run inversion and extract single models
model = INV.run()
res = model(0, f.nModel)  # first model part
vel = 1.0 / model(f.nModel, f.nModel*2)  # second model part
# %% show single models
pg.show(f.pd, res)
pg.show(f.pd, vel)
