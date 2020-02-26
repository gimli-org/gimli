#!/usr/bin/env python
# encoding: UTF-8

"""
"""

import pygimli as pg
import pybert as pb


class DCMultiElectrodeModellingC(pb.DCMultiElectrodeModelling):
    def __init__(self, mesh, data, verbose):
        super().__init__(mesh, data, verbose)
        self.setComplex(True)

        self._J = pg.matrix.BlockMatrix()
        self.setJacobian(self._J)

        self._C = pg.matrix.BlockMatrix()
        self.setConstraints(self._C)
        self.matrixHeap = []
        # super().createConstraints = self.createConstraints

    def createDefaultStartModel(self):
        """
        """
        res = pb.getComplexData(self.data())
        parCount = self.regionManager().parameterCount()
        re = pg.Vector(parCount, pg.mean(pg.math.real(res)))
        im = pg.Vector(parCount, -pg.mean(pg.math.imag(res)))
        return pg.cat(re, im)

    def createJacobian(self, model):
        print('=' * 100)

        if self.complex():
            modelRe = model[0:int(len(model)/2)]
            modelIm = model[int(len(model)/2):len(model)]
            modelC = pg.math.toComplex(modelRe, modelIm)
            print("Real", min(modelRe), max(modelRe))
            print("Imag", min(modelIm), max(modelIm))

            u = self.prepareJacobian_(modelC)

            if self._J.rows() == 0:
                #re(data)/re(mod) = im(data)/im(mod)
                # we need a local copy until we have a gimli internal reference counter FIXTHIS
                M1 = pg.Matrix()
                M2 = pg.Matrix()
                self.matrixHeap.append(M1)
                self.matrixHeap.append(M2)

                JRe = self._J.addMatrix(M1)
                JIm = self._J.addMatrix(M2)

                self._J.addMatrixEntry(JRe, 0, 0)
                self._J.addMatrixEntry(JIm, 0, len(modelRe), -1.0)
                self._J.addMatrixEntry(JIm, self.data().size(), 0, 1.0)
                self._J.addMatrixEntry(JRe, self.data().size(), len(modelRe))

            else:
                self._J.clean()


            k = pg.Vector(self.data()('k'))
            self.data().set('k', k*0.0 + 1.0)

            dMapResponse = pb.DataMap()
            dMapResponse.collect(self.electrodes(), self.solution())
            respRe = dMapResponse.data(self.data(), False, False)
            respIm = dMapResponse.data(self.data(), False, True)

            #CVector resp(toComplex(respRe, respIm));
            #RVector am(abs(resp) * dataContainer_->get("k"));
            #RVector ph(-phase(resp));

            print("respRe", pg.math.median(respRe), min(respRe), max(respRe))
            print("respIm", pg.math.median(respIm), min(respIm), max(respIm))

            JC = pg.matrix.CMatrix()
            self.createJacobian_(modelC, u, JC)
            for i in range(JC.rows()):
                #JC[i] *= 1.0/(modelC*modelC) * k[i]
                JC[i] /= (modelC * modelC) / k[i]

            self._J.mat(0).copy(pg.math.real(JC))
            self._J.mat(1).copy(pg.math.imag(JC))

            #self.createJacobian_(modelRe*0.0+1.0, pg.math.real(u), self._J.mat(1))
            #self.createJacobian_(modelRe*0.0+1.0, pg.math.imag(u), self._J.mat(2))
            #self.createJacobian_(modelRe*0.0+1.0, pg.math.imag(u), self._J.mat(3))


            sumsens0 = pg.Vector(self._J.mat(0).rows())
            sumsens1 = pg.Vector(self._J.mat(0).rows())
            sumsens2 = pg.Vector(self._J.mat(0).rows())

            for i in range(self._J.mat(0).rows()):
                #self._J.mat(0)[i] *= 1./modelRe / respRe[i]
                #self._J.mat(1)[i] *= 1./modelIm / respRe[i]

                #self._J.mat(2)[i] *= 1./modelRe / respIm[i]
                #self._J.mat(3)[i] *= 1./modelIm / respIm[i]

                #self._J.mat(0)[i] *= 1./(modelRe * modelRe) * k[i]
                #self._J.mat(1)[i] *= 1./(modelRe * modelIm) * k[i]

                #self._J.mat(2)[i] *= 1./(modelIm * modelRe) * k[i]
                #self._J.mat(3)[i] *= 1./(modelIm * modelIm) * k[i]

                sumsens0[i] = sum(self._J.mat(0)[i])
                sumsens1[i] = sum(self._J.mat(1)[i])
                sumsens2[i] = abs(sum(JC[i]))

            print(pg.math.median(sumsens0), min(sumsens0), max(sumsens0))
            print(pg.math.median(sumsens1), min(sumsens1), max(sumsens1))
            print(pg.math.median(sumsens2), min(sumsens2), max(sumsens2))

            self.data().set('k', k)

            self._J.recalcMatrixSize()
        else:
            # self.setVerbose(True)
            u = self.prepareJacobian_(model)

            #J = pg.Matrix()
            if self._J.rows() == 0:
                print('#' * 100)
                M1 = pg.Matrix()
                Jid = self._J.addMatrix(M1)
                self._J.addMatrixEntry(Jid, 0, 0)
            else:
                self._J.clean()

            self.createJacobian_(model, u, self._J.mat(0))
            self._J.recalcMatrixSize()

    def createConstraints(self):
        """
        """
        print ("createConstrains(self, model):")
        Ctmp = pg.matrix.SparseMapMatrix()
        self.matrixHeap.append(Ctmp)
        self.regionManager().fillConstraints(Ctmp)
        CiD = self._C.addMatrix(Ctmp)
        self._C.addMatrixEntry(CiD, 0, 0)
        self._C.addMatrixEntry(CiD, Ctmp.rows(), Ctmp.cols())
        self._C.recalcMatrixSize()


data = pb.DataContainerERT('wa24c.dat')
print(data)

mesh = pg.meshtools.createParaMesh2dGrid(data.sensorPositions())


fop = DCMultiElectrodeModellingC(mesh, data, verbose=True)
print(dir(fop))
print(fop.jacobian())
print(fop.jacobian().rows())

#fop = pb.DCMultiElectrodeModelling(mesh, data)

fop.regionManager().region(1).setBackground(True)
fop.createRefinedForwardMesh(refine=True, pRefine=False)

cData = pb.getComplexData(data)
mag = pg.abs(cData)
phi = -pg.phase(cData)

print(pg.norm(mag-data('rhoa')))
print(pg.norm(phi-data('ip')/1000))

inv = pg.Inversion(pg.cat(mag, phi),
                    fop,
                    verbose=True, dosave=True)


dataTrans = pg.trans.TransCumulative()
datRe = pg.trans.TransLog()
datIm = pg.trans.Trans()
dataTrans.add(datRe, data.size())
dataTrans.add(datIm, data.size())

modRe = pg.trans.TransLog()
modIm = pg.trans.TransLog()
modelTrans = pg.trans.TransCumulative()
modelTrans.add(modRe, fop.regionManager().parameterCount())
modelTrans.add(modIm, fop.regionManager().parameterCount())

inv.setTransData(dataTrans)
inv.setTransModel(modelTrans)
inv.setAbsoluteError(pg.cat(data("err")*mag, mag*phi*10.01))
inv.setLambda(5)
inv.setMaxIter(5)

#inv.setCWeight(pg.cat(pg.Vector(503, 1.0),
                      #pg.Vector(503, 100.0)))

model = inv.run()
jacRe = fop.jacobian()[0]
jacIm = fop.jacobian()[data.size()]

modelMesh = fop.regionManager().paraDomain()

pg.showLater(1)


fig = pg.plt.figure()
ax1 = fig.add_subplot(2, 2, 1)
ax2 = fig.add_subplot(2, 2, 2)
ax3 = fig.add_subplot(2, 2, 3)
ax4 = fig.add_subplot(2, 2, 4)
#pg.show(modelMesh, jacRe[0:len(jacRe)/2], colorBar=1, axes=ax1)
#ax1.set_title("jac real/real")
#pg.show(modelMesh, jacRe[len(jacRe)/2:len(jacRe)], colorBar=1, axes=ax2)
#ax2.set_title("jac real/imag")
#pg.show(modelMesh, jacIm[0:len(jacIm)/2], colorBar=1, axes=ax3)
#ax3.set_title("jac imag/real")
#pg.show(modelMesh, jacIm[len(jacIm)/2:len(jacIm)], colorBar=1, axes=ax4)
#ax4.set_title("jac imag/imag")


ax, cbar = pg.show(modelMesh, model[0:len(model)/2], colorBar=1)
ax.set_title("model real")
ax, cbar = pg.show(modelMesh, model[len(model)/2:len(model)], colorBar=1)
ax.set_title("model imag")


#
#pg.show(modelMesh, data.sensorPositions(), axes=ax, showLater=1)
#pg.show(modelMesh, axes=ax)
pg.showLater(0)