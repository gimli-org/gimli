#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Example for Petrophysical-Joint-Inversion."""

import numpy as np
import pygimli as pg

from pygimli import meshtools as mt
from pygimli.physics.petro import InvertPetro, InvertJointPetro
from pygimli.physics.petro import transFwdArchieS as ArchieTrans
from pygimli.physics.petro import transFwdWyllieS as WyllieTrans

def createSynthModel():
    """Return the modeling mesh, the porosity distribution and the
       parametric mesh for inversion.
    """
    ### Create the synthetic model
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

    ### Create the parametric mesh that only reflect the domain geometry
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


### Script starts here
axs = [None]*8

### Create synthetic model
mMesh, pMesh, saturation = createSynthModel()

### Create Petrophysical models
ertTrans = ArchieTrans(rFluid=20, phi=0.3)
res = ertTrans(saturation)

ttTrans = WyllieTrans(vm=4000, phi=0.3)
vel = 1./ttTrans(saturation)

### Simulate synthetic data with appropriate noise
sensors = mMesh.positions()[mMesh.findNodesIdxByMarker(-99)]

print("-Simulate ERT" + "-" * 50)
ERT = pg.physics.ERTManager(verbose=False)
ertScheme = pg.physics.ert.createERTData(sensors, schemeName='dd', closed=1)
ertData = ERT.simulate(mMesh, res, ertScheme, noiseLevel=0.01)

print("-Simulate Traveltime" + "-" * 50)
TT = pg.physics.Refraction()
ttScheme = pg.physics.traveltime.createRAData(sensors)
ttData = TT.simulate(mMesh, vel, ttScheme, noiseLevel=0.01, noiseAbs=4e-6)

## Classic inversions
print("-ERT" + "-" * 50)
resInv = ERT.invert(ertData, mesh=pMesh, zWeight=1, lam=20)
ERT.inv.echoStatus()

print("-TT" + "-" * 50)
velInv = TT.invert(ttData, mesh=pMesh, lam=100, useGradient=0, zWeight=1)
TT.inv.echoStatus()

### Petrophysical inversion
print("-ERT-Petro" + "-" * 50)
invERTPetro = InvertPetro(ERT, ertTrans)
satERT = invERTPetro.invert(ertData, mesh=pMesh, limits=[0., 1.], lam=10)
invERTPetro.inv.echoStatus()

print("-TT-Petro" + "-" * 50)
invTTPetro = InvertPetro(TT, ttTrans)
satTT = invTTPetro.invert(ttData, mesh=pMesh, limits=[0., 1.], lam=5)
invTTPetro.inv.echoStatus()

### Petrophysical joint inversion
print("-Joint-Petro" + "-" * 50)
invJointPetro = InvertJointPetro([ERT, TT], [ertTrans, ttTrans])
satJoint = invJointPetro.invert([ertData, ttData], mesh=pMesh, limits=[0., 1.], lam=5)
invJointPetro.inv.echoStatus()

### Show results
ERT.showData(ertData)
TT.showVA(ttData)

showModel(axs[0], saturation, mMesh, showMesh=1,
          label=r'Saturation (${\tt petro}$)', savefig='petro')
showModel(axs[1], res, mMesh, petro=0, cMin=250, cMax=2500, showMesh=1,
          label=r'Resistivity (${\tt res}$) in $\Omega$m', savefig='resistivity')
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

# just hold the figure windows open if run outside from spyder, ipython or similar
pg.wait()
