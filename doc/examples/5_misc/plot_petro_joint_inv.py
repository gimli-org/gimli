#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Petrophysical joint inversion
-----------------------------

Joint inversion of different geophysical techniques helps to improve both
resolution and interpretability of the resulting images. Different data sets can
be directly coupled, if there is a link to an underlying target parameter.
In this example, ERT and traveltime data are inverted for water saturation. For
details see section 3.3 of the pyGIMLi paper (https://cg17.pygimli.org).
"""
# sphinx_gallery_thumbnail_number = 5

import numpy as np
import pygimli as pg

from pygimli import meshtools as mt
from pygimli.physics.petro import transFwdArchieS as ArchieTrans
from pygimli.physics.petro import transFwdWyllieS as WyllieTrans
from pygimli.frameworks import PetroInversionManager, JointPetroInversionManager


################################################################################
# We start with defining two helper functions.

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

    # Create the parametric mesh that only reflect the domain geometry
    world = mt.createCircle(boundaryMarker=-1, segments=32, area=0.0051)
    paraMesh = pg.meshtools.createMesh(world, q=34.0, smooth=[1, 10])
    paraMesh.scale(1.0/5.0)

    return mesh, paraMesh, petro


def showModel(ax, model, mesh, petro=1, cMin=None, cMax=None, label=None,
              cMap=None, showMesh=False):
    """Utility function to show and save models for the CG paper."""
    if cMin is None:
        cMin = 0.3
    if cMax is None:
        cMax = 0.9

    if cMap is None:
        cMap = 'viridis'
    if petro:
        ax, _ = pg.show(mesh, model, label=label,
                        logScale=False, cMin=cMin, cMax=cMax, cMap=cMap, ax=ax)
    else:
        ax, _ = pg.show(mesh, model, label=label,
                        logScale=True, cMin=cMin, cMax=cMax, cMap=cMap, ax=ax)

    ticks = [-.2, -.1, 0, .1, .2]
    ax.xaxis.set_ticks(ticks)
    ax.yaxis.set_ticks(ticks)

    pg.viewer.mpl.drawSensors(ax, ertData.sensorPositions(), diam=0.005)

    # despine(ax=ax, offset=5, trim=True)
    if showMesh:
        pg.viewer.mpl.drawSelectedMeshBoundaries(ax, mesh.boundaries(),
                                                 linewidth=0.3, color="0.2")
    return ax

################################################################################
# Create synthetic model
# ......................
mMesh, pMesh, saturation = createSynthModel()

################################################################################
# Create Petrophysical models
ertTrans = ArchieTrans(rFluid=20, phi=0.3)
res = ertTrans(saturation)

ttTrans = WyllieTrans(vm=4000, phi=0.3)
vel = 1./ttTrans(saturation)

sensors = mMesh.positions()[mMesh.findNodesIdxByMarker(-99)]

################################################################################
# Forward simulation
# ..................
# To create synthetic data sets, we assume 16 equally-spaced sensors on the
# circumferential boundary of the mesh. For the ERT modelling we build a
# complete dipole-dipole array. For the ultrasonic tomography we simulate the
# travel time for every possible sensor pair.
pg.info("Simulate ERT")
ERT = pg.physics.ert.ERTManager(verbose=False, sr=False)
ertScheme = pg.physics.ert.createERTData(sensors, schemeName='dd', closed=1)
ertData = ERT.simulate(mMesh, scheme=ertScheme, res=res, noiseLevel=0.01)

pg.info("Simulate Traveltime")
TT = pg.physics.traveltime.TravelTimeManager(verbose=False)
ttScheme = pg.physics.traveltime.createRAData(sensors)
ttData = TT.simulate(mMesh, scheme=ttScheme, vel=vel,
                     noiseLevel=0.01, noiseAbs=4e-6)

################################################################################
# Conventional inversion
pg.info("ERT Inversion")
resInv = ERT.invert(ertData, mesh=pMesh, zWeight=1, lam=20, verbose=False)
ERT.inv.echoStatus()

pg.info("Traveltime Inversion")
velInv = TT.invert(ttData, mesh=pMesh, lam=100, useGradient=0, zWeight=1.0)
TT.inv.echoStatus()

################################################################################
# Petrophysical inversion (individually)
pg.info("ERT Petrogeophysical Inversion")
ERTPetro = PetroInversionManager(petro=ertTrans, mgr=ERT)
satERT = ERTPetro.invert(ertData, mesh=pMesh, limits=[0., 1.], lam=10,
                         verbose=False)
ERTPetro.inv.echoStatus()

pg.info("TT Petrogeophysical Inversion")
TTPetro = PetroInversionManager(petro=ttTrans, mgr=TT)
satTT = TTPetro.invert(ttData, mesh=pMesh, limits=[0., 1.], lam=5)
TTPetro.inv.echoStatus()

################################################################################
# Petrophysical joint inversion
pg.info("Petrophysical Joint-Inversion TT-ERT")
JointPetro = JointPetroInversionManager(petros=[ertTrans, ttTrans],
                                        mgrs=[ERT, TT])
satJoint = JointPetro.invert([ertData, ttData], mesh=pMesh,
                             limits=[0., 1.], lam=5, verbose=False)
JointPetro.inv.echoStatus()

################################################################################
# Show results
ERT.showData(ertData)
TT.showData(ttData)

axs = [None]*8

showModel(axs[0], saturation, mMesh, showMesh=True,
          label=r'Saturation (${\tt petro}$)')
showModel(axs[1], res, mMesh, petro=0, cMin=250, cMax=2500, showMesh=1,
          label=pg.unit('res'), cMap=pg.cmap('res'))
showModel(axs[5], vel, mMesh, petro=0, cMin=1000, cMax=2500, showMesh=1,
          label=pg.unit('vel'), cMap=pg.cmap('vel'))
showModel(axs[2], resInv, pMesh, 0, cMin=250, cMax=2500,
          label=pg.unit('res'), cMap=pg.cmap('res'))
showModel(axs[6], velInv, pMesh, 0, cMin=1000, cMax=2500,
          label=pg.unit('vel'), cMap=pg.cmap('vel'))
showModel(axs[3], satERT, pMesh,
          label=r'Saturation (${\tt satERT}$)')
showModel(axs[7], satTT, pMesh,
          label=r'Saturation (${\tt satTT}$)')
showModel(axs[4], satJoint, pMesh,
          label=r'Saturation (${\tt satJoint}$)')
