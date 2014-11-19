#!/usr/bin/env python
# encoding: UTF-8

"""
"""

import pygimli as pg
import pybert as pb


class BlockMatrix(pg.MatrixBase):
    def __init__(self):
        super().__init__()
        self._matrieces = []
        self._indieces = []
        self._rows = 0
        self._cols = 0

    def rows(self):
        return self._rows

    def cols(self):
        return self._cols

    def row(self, idx):
        ret = pg.RVector(self._cols)
        for rowStart, colStart, matrixID, scale in self._indieces:
            mat = self._matrieces[matrixID]
            # print(mat[0], mat.cols(), mat.rows(), rowStart)
            if idx >= rowStart and idx < rowStart + mat.rows():
                # print(mat[idx - rowStart], mat.cols(), mat.rows())
                ret.setVal(mat.row(idx - rowStart) * scale,
                           colStart, colStart + mat.cols())
        return ret

    def col(self, idx):
        ret = pg.RVector(self._rows)
        for rowStart, colStart, matrixID, scale in self._indieces:
            mat = self._matrieces[matrixID]
            # print(mat[0], mat.cols(), mat.rows(), rowStart)
            if idx >= colStart and idx < colStart + mat.cols():
                ret.setVal(mat.col(idx - colStart) * scale,
                           rowStart, rowStart + mat.rows())
        return ret
    
    def __getitem__(self, idx):
        if type(idx) == int:
            return self.row(idx)
        raise

    def clear(self):
        self._matrieces = []
        self._indieces = []
    def clean(self):
        for m in self._matrieces:
            m.clean()
    
    def mat(self, idx):
        return self._matrieces[idx]

    def addMatrix(self, matrix):
        self._matrieces.append(matrix)
        return len(self._matrieces) - 1

    def addMatrixEntry(self, matrixID, rowStart, colStart, scale=1.0):
        if matrixID > len(self._matrieces):
            raise
        self._indieces.append([rowStart, colStart, matrixID, scale])
        self._rows = max(self._rows,
                         rowStart + self._matrieces[matrixID].rows())
        self._cols = max(self._cols,
                         colStart + self._matrieces[matrixID].cols())

    def recalcMatrixSize(self):
        for rowStart, colStart, matrixID, scale in self._indieces:
            self._rows = max(self._rows,
                             rowStart + self._matrieces[matrixID].rows())
            self._cols = max(self._cols,
                             colStart + self._matrieces[matrixID].cols())

    def mult(self, b):
        if len(b) != self.cols():
            raise Exception("matrix mult dimension missmatch" +
                            str(self.cols()) + " " + str(len(b)))

        ret = pg.RVector(self._rows)
        
        for rowStart, colStart, matrixID, scale in self._indieces:
            mat = self._matrieces[matrixID]
            
            ret.addVal(mat.mult(b.getVal(colStart, colStart + mat.cols())) * scale,
                       rowStart, rowStart + mat.rows())
        
        return ret

    def transMult(self, b):
        if len(b) != self.rows():
            raise Exception("matrix mult dimension missmatch" +
                            str(self.rows()) + " " + str(len(b)))

        ret = pg.RVector(self._cols)
        
        for rowStart, colStart, matrixID, scale in self._indieces:
            mat = self._matrieces[matrixID]
            
            ret.addVal(mat.transMult(b.getVal(rowStart, rowStart + mat.rows())) * scale,
                       colStart, colStart + mat.cols())
        return ret

    def save(self, name):
        print("won't save matrix")
        


class DCMultiElectrodeModellingC(pb.DCMultiElectrodeModelling):
    def __init__(self, mesh, data, verbose):
        super().__init__(mesh, data, verbose)
        self.setComplex(True)
        
        self._J = BlockMatrix()
        self.setJacobian(self._J)
    
        self._C = pg.RBlockMatrix()
        self.setConstraints(self._C)
        #super().createConstraints = self.createConstraints
    
    def createDefaultStartModel(self):
        """
        """
        res = pb.getComplexData(self.data())
        re = pg.RVector(self.regionManager().parameterCount(), pg.mean(pg.real(res)))
        im = pg.RVector(self.regionManager().parameterCount(), pg.mean(pg.imag(res)))
        return pg.cat(re, im)
        
    def createJacobian(self, model):
        print('=' * 100)
                
        if self.complex():
            modelRe = model[0:int(len(model)/2)]
            modelIm = model[int(len(model)/2):len(model)]
            modelC = pg.toComplex(modelRe, modelIm)
            print("Real", min(modelRe), max(modelRe))
            print("Imag", min(modelIm), max(modelIm))
            
            u = self.prepareJacobian_(modelC)
            
            if self._J.rows() == 0:
                #re(data)/re(mod) = im(data)/im(mod)
                # we need a local copy until we have a gimli internal reference counter FIXTHIS
                M1 = pg.RMatrix()
                M2 = pg.RMatrix()
                M3 = pg.RMatrix()
                M4 = pg.RMatrix()
                JReReID = self._J.addMatrix(M1)
                JImReID = self._J.addMatrix(M2)
                JReImID = self._J.addMatrix(M3)
                JImImID = self._J.addMatrix(M4)
                
                self._J.addMatrixEntry(JReReID, 0, 0)
                self._J.addMatrixEntry(JImReID, 0, len(modelRe))
                self._J.addMatrixEntry(JReImID, self.data().size(), 0, -1)
                self._J.addMatrixEntry(JImImID, self.data().size(), len(modelRe))
                
            else:
                self._J.clean()
                
                
            k = pg.RVector(self.data()('k'))
            self.data().set('k', k*0.0 + 1.0)
            
            dMapResponse = pb.DataMap()
            dMapResponse.collect(self.electrodes(), self.solution())
            respRe = dMapResponse.data(self.data(), False, False)
            respIm = dMapResponse.data(self.data(), False, True)
                
            #CVector resp(toComplex(respRe, respIm));
            #RVector am(abs(resp) * dataContainer_->get("k"));
            #RVector ph(-phase(resp));
                
            print("respRe", pg.median(respRe), min(respRe), max(respRe))
            print("respIm", pg.median(respIm), min(respIm), max(respIm))                
            
            self.createJacobian_(modelRe*0.0+1.0, pg.real(u), self._J.mat(0))
            self.createJacobian_(modelRe*0.0+1.0, pg.real(u), self._J.mat(1))
            
            self.createJacobian_(modelRe*0.0+1.0, pg.imag(u), self._J.mat(2))
            self.createJacobian_(modelRe*0.0+1.0, pg.imag(u), self._J.mat(3))
                       
            
            sumsens0 = pg.RVector(self._J.mat(0).rows())
            sumsens1 = pg.RVector(self._J.mat(0).rows())
            sumsens2 = pg.RVector(self._J.mat(0).rows())
            sumsens3 = pg.RVector(self._J.mat(0).rows())
            for i in range(self._J.mat(0).rows()):
                #self._J.mat(0)[i] *= 1./modelRe / respRe[i]
                #self._J.mat(1)[i] *= 1./modelIm / respRe[i]
                
                #self._J.mat(2)[i] *= 1./modelRe / respIm[i]
                #self._J.mat(3)[i] *= 1./modelIm / respIm[i]
                
                self._J.mat(0)[i] *= 1./(modelRe * modelRe) * k[i]
                self._J.mat(1)[i] *= 1./(modelRe * modelIm) * k[i]
                
                self._J.mat(2)[i] *= 1./(modelIm * modelRe) * k[i]
                self._J.mat(3)[i] *= 1./(modelIm * modelIm) * k[i]
                
                sumsens0[i] = sum(self._J.mat(0)[i])
                sumsens1[i] = sum(self._J.mat(1)[i])
                sumsens2[i] = sum(self._J.mat(2)[i])
                sumsens3[i] = sum(self._J.mat(3)[i])
                
            print(pg.median(sumsens0), min(sumsens0), max(sumsens0))
            print(pg.median(sumsens1), min(sumsens1), max(sumsens1))
            print(pg.median(sumsens2), min(sumsens2), max(sumsens2))
            print(pg.median(sumsens3), min(sumsens3), max(sumsens3))
            
            #for i in range(self._J.mat(0).rows()):
                #print(pg.norm(self._J.mat(0)[i] - self._J.mat(1)[i]))
                #print(pg.norm(self._J.mat(2)[i] - self._J.mat(3)[i]))
                
            self.data().set('k', k)
            
            self._J.recalcMatrixSize()
        else:
            # self.setVerbose(True)
            u = self.prepareJacobian_(model)
            
            #J = pg.RMatrix()
            if self._J.rows() == 0:
                print('#' * 100)
                M1 = pg.RMatrix()
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
        Ctmp = pg.RSparseMapMatrix()
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

inv = pg.RInversion(pg.cat(mag, phi),
                    fop, verbose=True, dosave=True)

dataTrans = pg.RTransCumulative()
datRe = pg.RTransLog()
datIm = pg.RTrans()
dataTrans.add(datRe, data.size())
dataTrans.add(datIm, data.size())

modRe = pg.RTransLog()
modIm = pg.RTrans()
modelTrans = pg.RTransCumulative()
modelTrans.add(modRe, fop.regionManager().parameterCount())
modelTrans.add(modIm, fop.regionManager().parameterCount())

inv.setTransData(dataTrans)
inv.setTransModel(modelTrans)
inv.setError(pg.cat(data("err"), data("err")))
inv.setLambda(5)
inv.setMaxIter(5)

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
pg.show(modelMesh, jacRe[0:len(jacRe)/2], colorBar=1, axes=ax1)
ax1.set_title("jac real/real")
pg.show(modelMesh, jacRe[len(jacRe)/2:len(jacRe)], colorBar=1, axes=ax2)
ax2.set_title("jac real/imag")
pg.show(modelMesh, jacIm[0:len(jacIm)/2], colorBar=1, axes=ax3)
ax3.set_title("jac imag/real")
pg.show(modelMesh, jacIm[len(jacIm)/2:len(jacIm)], colorBar=1, axes=ax4)
ax4.set_title("jac imag/imag")


ax, cbar = pg.show(modelMesh, model[0:len(model)/2], colorBar=1)
ax.set_title("model real")
ax, cbar = pg.show(modelMesh, model[len(model)/2:len(model)], colorBar=1)
ax.set_title("model imag")


#
#pg.show(modelMesh, data.sensorPositions(), axes=ax, showLater=1)
#pg.show(modelMesh, axes=ax)
pg.showLater(0)