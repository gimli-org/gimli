#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Example for Petrophysical-Joint-Inversion.

The presented classes are basic versions for demonstration purposes and
without any further syntactic sugar or any input or probability checks.

More advanced version can be find in the pygimli library itself.
"""

import numpy as np
import pygimli as pg

from pygimli import meshtools as mt

from pygimli.frameworks import MeshInversion
from pygimli.frameworks import Modelling

from pygimli.physics.petro import transFwdArchieS as ArchieTrans
from pygimli.physics.petro import transFwdWyllieS as WyllieTrans


def createSynthModel():
    """Return the modelling mesh, the porosity distribution and the
       parametric mesh for inversion.
    """
    # Create the synthetic model
    world = mt.createCircle(boundaryMarker=-1, segments=64)
    tri = mt.createPolygon([[-0.8, -0], [-0.5, -0.7], [0.7, 0.5]],
                           isClosed=True, area=0.0015)
    c1 = mt.createCircle(radius=0.2, pos=[-0.2, 0.5], segments=32,
                         area=0.0025, marker=3)
    c2 = mt.createCircle(radius=0.2, pos=[0.32, -0.3], segments=32,
                         area=0.0025, marker=3)

    poly = mt.mergePLC([world, tri, c1, c2])

    poly.addRegionMarker([0.0, 0, 0], 1, area=0.0015)
    poly.addRegionMarker([-0.9, 0, 0], 2, area=0.0015)

    c = mt.createCircle(radius=0.99, segments=16, start=np.pi, end=np.pi*3)
    [poly.createNode(p.pos(), -99) for p in c.nodes()]
    mesh = pg.meshtools.createMesh(poly, q=34.4, smooth=[1, 10])
    mesh.scale(1.0/5.0)
    mesh.rotate([0., 0., 3.1415/3])
    mesh.rotate([0., 0., 3.1415])

    petro = pg.solver.parseArgToArray([[1, 0.9], [2, 0.6], [3, 0.3]],
                                      mesh.cellCount(), mesh)

    # Create the parametric mesh that only reflects the domain geometry
    world = mt.createCircle(boundaryMarker=-1, segments=32, area=0.0051)
    paraMesh = pg.meshtools.createMesh(world, q=34.0, smooth=[1, 10])
    paraMesh.scale(1.0/5.0)

    return mesh, paraMesh, petro


def showModel(ax, model, mesh, petro=1, cMin=None, cMax=None, label=None,
              savefig=None, showMesh=False):
    """Utility function to show and save models for the CG paper."""
    if cMin is None:
        cMin = 0.3
    if cMax is None:
        cMax = 1.0

    if petro:
        ax, _ = pg.show(mesh, model, label=label,
                        cMin=cMin, cMax=cMax, logScale=0, ax=ax,
                        cMap='viridis', hold=1)
    else:
        ax, _ = pg.show(mesh, model, label=label,
                        logScale=1, ax=ax, cMin=cMin, cMax=cMax, hold=1)

    ticks = [-.2, -.1, 0, .1, .2]
    ax.xaxis.set_ticks(ticks)
    ax.yaxis.set_ticks(ticks)

    pg.mplviewer.drawSensors(ax, ertData.sensorPositions(), diam=0.005)

    # despine(ax=ax, offset=5, trim=True)
    if showMesh:
        pg.mplviewer.drawSelectedMeshBoundaries(ax, mesh.boundaries(),
                                                linewidth=0.3, color="0.2")

    if savefig:
        pg.mplviewer.saveAxes(ax, savefig, adjust=False)
    return ax


class PetroModelling(Modelling):
    """Combine relation m(p) with modelling f(p). """
    def __init__(self, fop, petro):
        """Initialize with instance of forward operator
        and transformation, create a MultRightMatrix Jacobian."""
        Modelling.__init__(self)
        self.fop = fop
        self.petro = petro
        self.jac = pg.MultRightMatrix(self.fop.jacobian())
        self.setJacobian(self.jac)

    def response(self, model):
        """Transform and compute response f(p(m))."""
        tModel = self.petro.fwd(model)
        return self.fop.response(tModel)

    def createJacobian(self, model):
        """Jacobian with inner derivative of the trans."""
        self.fop.createJacobian(self.petro.fwd(model))
        self.jac.r = self.petro.deriv(model)


class JointModelling(Modelling):
    """Cumulative (joint) forward operator."""
    def __init__(self, fopList):
        """Initialize with lists of forward operators"""
        Modelling.__init__(self)
        self.fops = fopList
        self.jac = pg.matrix.BlockMatrix()

    def response(self, model):
        """Concatenate responses for all fops."""
        resp = []
        for f in self.fops:
            resp.extend(f.response(model))
        return resp

    def createJacobian(self, model):
        """Fill the individual Jacobian matrices."""
        self.initJacobian()
        for f in self.fops:
            f.createJacobian(model)

    def setData(self, data):
        nData = 0
        for i, fi in enumerate(self.fops):
            fi.setData(data[i])
            self.jac.addMatrix(fi.jacobian(), nData, 0)
            nData += data[i].size()  # update total vector length
        self.setJacobian(self.jac)

    def setMesh(self, mesh):
        for fi in self.fops:
            fi.setMesh(mesh)
        self.setRegionManager(self.fops[0].regionManagerRef())


class PetroInversion(MeshInversion):
    def __init__(self, mgr, petro):
        MeshInversion.__init__(self)
        self.mgr = mgr
        self.fop = PetroModelling(mgr.createFOP(), petro)
        self.tM = mgr.tM
        self.tD = mgr.tD
        self.inv.setTransData(self.tD)
        self.setForwardOperator(self.fop)

    def setData(self, data):
        self.fop.setData(data)
        self.dataVals = self.mgr.dataVals(data)
        self.dataErrs = self.mgr.relErrorVals(data)

    def invert(self, data, **kwargs):
        limits = kwargs.pop('limits', [0., 1.])
        self.tM.setLowerBound(limits[0])
        self.tM.setUpperBound(limits[1])
        self.inv.setTransModel(self.tM)
        kwargs['startModel'] = (limits[1]-limits[0]) / 2.
        return MeshInversion.invert(self, data, **kwargs)


class JointPetroInversion(MeshInversion):
    def __init__(self, mgrs, petros):
        """Initialize with lists of managers and transformations"""
        MeshInversion.__init__(self)
        self.mgrs = mgrs

        self.fops = [PetroModelling(m.createFOP(), p)
                     for m, p in zip(mgrs, petros)]

        self.tM = mgrs[0].tM
        self.tD = pg.trans.TransCumulative()
        self.fop = JointModelling(self.fops)
        self.setForwardOperator(self.fop)

    def setData(self, data):
        self.fop.setData(data)
        self.dataVals = pg.Vector(0)
        self.dataErrs = pg.Vector(0)

        for i, mgr in enumerate(self.mgrs):
            self.tD.add(mgr.tD, data[i].size())
            self.dataVals = pg.cat(self.dataVals, mgr.dataVals(data[i]))
            self.dataErrs = pg.cat(self.dataErrs, mgr.relErrorVals(data[i]))
        self.inv.setTransData(self.tD)

    def invert(self, data, **kwargs):
        limits = kwargs.pop('limits', [0., 1.])
        self.tM.setLowerBound(limits[0])
        self.tM.setUpperBound(limits[1])
        self.inv.setTransModel(self.tM)
        kwargs['startModel'] = (limits[1]-limits[0])/2.
        return MeshInversion.invert(self, data, **kwargs)

# Script starts here #
# Create synthetic model
mMesh, pMesh, saturation = createSynthModel()

# Create Petrophysical models
ertTrans = ArchieTrans(rFluid=20, phi=0.3)
res = ertTrans(saturation)

ttTrans = WyllieTrans(vm=4000, phi=0.3)
vel = 1./ttTrans(saturation)

# Simulate synthetic data with appropriate noise
sensors = mMesh.positions()[mMesh.findNodesIdxByMarker(-99)]

print("-Simulate ERT" + "-" * 50)
ERT = pg.physics.ERTManager(verbose=False)
ertScheme = pg.physics.ert.createERTData(sensors, schemeName='dd', closed=1)
ertData = ERT.simulate(mMesh, res, ertScheme, noiseLevel=0.01)

print("-Simulate Traveltime" + "-" * 50)
TT = pg.physics.Refraction(verbose=False)
ttScheme = pg.physics.traveltime.createRAData(sensors)
ttData = TT.simulate(mMesh, vel, ttScheme, noiseLevel=0.01, noiseAbs=4e-6)

# Classic inversions
print("-ERT" + "-" * 50)
resInv = ERT.invert(ertData, mesh=pMesh, zWeight=1, lam=20)
ERT.inv.echoStatus()

print("-TT" + "-" * 50)
velInv = TT.invert(ttData, mesh=pMesh, lam=100, useGradient=0, zWeight=1)
TT.inv.echoStatus()

print("-ERT-Petro" + "-" * 50)
invERTPetro = PetroInversion(ERT, ertTrans)
satERT = invERTPetro.invert(ertData, mesh=pMesh, limits=[0., 1.], lam=10)
invERTPetro.inv.echoStatus()

print("-TT-Petro" + "-" * 50)
invTTPetro = PetroInversion(TT, ttTrans)
satTT = invTTPetro.invert(ttData, mesh=pMesh, limits=[0., 1.], lam=5)
invTTPetro.inv.echoStatus()


# Petrophysical joint inversion
print("-Joint-Petro" + "-" * 50)
invJointPetro = JointPetroInversion([ERT, TT], [ertTrans, ttTrans])
satJoint = invJointPetro.invert([ertData, ttData], mesh=pMesh,
                                limits=[0., 1.], lam=5)
invJointPetro.inv.echoStatus()

# Show results
ERT.showData(ertData)
TT.showVA(ttData)

axs = [None]*8
showModel(axs[0], saturation, mMesh, showMesh=1,
          label=r'Saturation (${\tt petro}$)', savefig='petro')
showModel(axs[1], res, mMesh, petro=0, cMin=250, cMax=2500, showMesh=1,
          label=r'Resistivity (${\tt res}$) in $\Omega$m',
          savefig='resistivity')
showModel(axs[5], vel, mMesh, petro=0, cMin=1000, cMax=2500, showMesh=1,
          label=r'Velocity (${\tt vel}$) in m$/$s', savefig='velocity')
showModel(axs[2], resInv, pMesh, 0, cMin=250, cMax=2500,
          label=r'Resistivity (${\tt resInv}$) in $\Omega$m', savefig='invERT')
showModel(axs[6], velInv, pMesh, 0, cMin=1000, cMax=2500,
          label=r'Velocity (${\tt velInv}$) in m$/$s', savefig='invTT')
showModel(axs[3], satERT, pMesh,
          label=r'Saturation (${\tt satERT}$)', savefig='invERTPetro')
showModel(axs[7], satTT, pMesh,
          label=r'Saturation (${\tt satTT}$)', savefig='invTTPetro')
showModel(axs[4], satJoint, pMesh,
          label=r'Saturation (${\tt satJoint}$)', savefig='invJointPetro')

# just hold figure windows open if run outside from spyder, ipython or similar
pg.wait()
