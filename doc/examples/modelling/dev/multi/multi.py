#!/usr/bin/env python

"""
Test multi
"""

import os
import sys
import time

import numpy as np

import matplotlib
#matplotlib.use('pdf')
import matplotlib.animation as animation
import matplotlib.pyplot as plt

import pygimli as pg
import pygimli.polytools as pt
from pygimli.meshtools import createMesh

import fluidFlow
import seismics
import gravimetry
import ert

def savefig(mesh, plc, data, label, out=None, showMesh=0):
    if savefigs_:
        ax, _ = pg.show(mesh, data, hold=1, colorBar=1, label=label)
        if showMesh:
            pg.show(mesh, axes=ax, hold=1)
        if out:
            pg.show(plc, axes=ax, savefig=out + '.pdf')
            try:
                print("trying pdf2pdfS ... ")
                os.system('pdf2pdfS ' + out + '.pdf')
            except:
                pass
        else:
            pg.show(plc, axes=ax)

def createAnimation(fig, animate, nFrames, dpi, out):
    anim = animation.FuncAnimation(fig, animate,
                                   frames=nFrames,
                                   interval=0.001, repeat=False)
    anim.save(out + ".mp4", writer=None, fps=20, dpi=dpi, codec=None,
              bitrate=24*1024, extra_args=None, metadata=None,
              extra_anim=None, savefig_kwargs=None)
    try:
        print("create frames ... ")
        os.system('mkdir -p anim-' + out)
        os.system('ffmpeg -i ' + out + '.mp4 anim-' + out + '/movie%d.jpg')
    except:
        pass


def saveani(mesh, plc, data, label, out,
            cMin=None, cMax=None, logScale=False, cmap=None):
    """
    """
    dpi=92
    scale=1
    fig = plt.figure(facecolor='white',
                     figsize=(scale*800/dpi, scale*490/dpi), dpi=dpi)  
    ax = fig.add_subplot(1,1,1)
        
    gci = pg.mplviewer.drawModel(ax, mesh, data=data[0],
                                 cMin=cMin, cMax=cMax, cmap=cmap,
                                 logScale=logScale)
    
    cbar = pg.mplviewer.createColorbar(gci, label=label)
    pg.show(plc, axes=ax)
    
    plt.tight_layout()
    plt.pause(0.001)
    
    def animate(i):
        print(out + ": Frame:", i, "/", len(data))
        pg.mplviewer.setMappableData(gci, pg.abs(data[i]), 
                                     cMin=cMin, cMax=cMax,
                                     logScale=logScale)
        #plt.pause(0.001)
    createAnimation(fig, animate, int(len(data)), dpi, out)
    
    
   
    
def createTestWorld1(maxArea=0.2, verbose=0):
    
    #  ___________________8___________________
    #  |                                     |
    #  1                                     7
    #  |__________________9__________________|
    #  |                                     |
    #  2            |     4     |            6
    #  |__________________10_________________|
    #  |                                     |
    #  3                                     5 
    #  |__________________4__________________|

    layer1 = pt.createRectangle(start=[-20, 0], end=[20, -2], 
                                 marker=1, boundaryMarker=[1, 9, 7, 8])
    layer2 = pt.createRectangle(start=[-20, -2], end=[20, -8], 
                                 marker=2, boundaryMarker=[2, 10, 6, 9])
    layer3 = pt.createRectangle(start=[-20, -8], end=[20, -15], 
                                 marker=3, boundaryMarker=[3, 4, 5, 10])
    block = pt.createRectangle(start=[-6, -3.5], end=[6, -6.0], marker=4,
                               boundaryMarker=10)
    
    plc = pt.mergePLC([layer1, layer2, layer3, block])
    
    mesh = createMesh(plc, quality=33, area=maxArea,
                      smooth=[1,10], verbose=verbose)
   
    return mesh, plc
    
swatch = pg.Stopwatch(True)

mesh, plc = createTestWorld1(maxArea=0.2, verbose=0)

poro = pg.solver.parseMapToCellArray([[1, 0.3], 
                                      [2, 0.4], 
                                      [3, 0.2],
                                      [4, 0.115]], mesh)
velBoundary=[[2, [1.0, 0.0]],
             [6, [1.0, 0.0]],
             [3, [1./2., 0.0]],
             [8, [0.0, 0.0]],
             [5, [1./2, 0.0]],
             [4, [1./2, 0.0]]
             ]
preBoundary=[[1, 0.0],
             [7, 0.0]]

perm = fluidFlow.permeabiltyEngelhardtPitter(poro, mesh=mesh)    
hydr = fluidFlow.hydraulicConductivity(perm, mesh=mesh)
vP   = seismics.velocityVp(poro, mesh=mesh)
dens = gravimetry.density(poro, densMatrix=2510, densFluid=1000, mesh=mesh)

savefigs_ = 0
if savefigs_:
    savefig(mesh, plc, mesh.cellMarker(), 'Marker', 'marker')
    savefig(mesh, plc, mesh.cellMarker(), 'Marker', 'mesh', 1)
    savefig(mesh, plc, poro, 'Porosity', 'poro')
    savefig(mesh, plc, perm, 'Permeabilty [m$^2$]', 'perm')
    savefig(mesh, plc, hydr, 'Hydraulic conductivity [m$/$s]', 'hydr')
    savefig(mesh, plc, vP, 'Seismic velocity [m$/$s]', 'vP')
    savefig(mesh, plc, dens, 'Density [kg$/$m$^3$]', 'dens')

grav = 0
if grav:
    Grav, densBlock = gravimetry.calcInvBlock(mesh, dens, 'gravInv.pdf')
    savefig(mesh, plc, densBlock, 'Delta Density [kg$/$m$^3$]', 'ddens')
    savefig(Grav.fop.regionManager().paraDomain(), plc,
            Grav.inv.model(), 'Delta Density [kg$/$m$^3$]', 'densGrav')

seis = 0
if seis:
    seismics.calcSeismics(mesh, vP)

visc = 1./perm
vel = fluidFlow.calcStokesVelocity(mesh, visc, velBoundary, preBoundary)
print("vel:", swatch.duration(True))

showVel = 0
if showVel:
    axVel, _ = pg.show(mesh, data=pg.cellDataToPointData(mesh,
                              np.sqrt(vel[:,0]*vel[:,0] + vel[:,1]*vel[:,1])),
                       logScale=0, colorBar=1,
                       label='|Velocity| in m/s', hold=1)

    pg.viewer.showBoundaryNorm(mesh, velBoundary, color='red', axes=axVel)
    pg.show(plc, axes=axVel)
    
    meshC, _ = createTestWorld1(maxArea=1, verbose=0)
    pg.show(mesh, data=vel, axes=axVel, coarseMesh=meshC, savefig='velocity.pdf')
    

################## TIMELAPSE Starts here                  ######################
timeSteps = np.linspace(0, 50, 20) 

################## TIMELAPSE Concentration                ######################
conc = fluidFlow.calcConcentration(mesh, vel, timeSteps, 
                                   injectPos=[-18., -6], scale=1, peclet=50)

print("conc:", swatch.duration(True))


################## TIMELAPSE ERT                          ######################
rhoaMesh, _ = createTestWorld1(maxArea=0.1, verbose=0)

rhoaMesh, resis, ertData, rhoa = ert.calcApparentResistivities(mesh, 
                                                                 rhoaMesh, poro,
                                          rhoBrine=1./(1./20. + abs(1.*conc)))
ertMesh, ertMod, ertRat = ert.calcERT(ertData, rhoa)

print("ert:", swatch.duration(True))
saveanim = 0:
if saveanim:
    saveani(mesh, plc, conc, 'Concentration [kg/m$^3$]', 'conc', cMin=0, cMax=0.1)
    saveani(rhoaMesh, plc, resis, 'Resistivity [$\Omega$ m]', 'resis', logScale=1)
    saveani(ertMesh, plc, ertMod, 'Resistivity [$\Omega$ m]', 'ert', logScale=1, cMin=100, cMax=1000)
    saveani(ertMesh, plc, ertRat, 'Ratio', 'ratio', logScale=1, cMin=0.5, cMax=2, cmap='bwr')


# ** make Apparent resistivities animations (data) **
Ert = ert.ERT()
fig, ax = plt.subplots()
gci = Ert.show(ertData, values=rhoa[0], axes=ax, scheme='DipoleDipole', cMin=150, cMax=350)

def animate(i):
    print("rhoa: Frame:", i, "/", len(rhoa))
    ax.clear()
    Ert.show(ertData, values=rhoa[i], axes=ax, scheme='DipoleDipole', cMin=150, cMax=350, 
             colorBar=0)
    plt.tight_layout()
createAnimation(fig, animate, int(len(rhoa)), 92, 'rhoa')   

fig, ax = plt.subplots()
gci = Ert.show(ertData, values=rhoa[0], axes=ax, scheme='DipoleDipole', cMin=0.5, cMax=2)
def animate(i):
    print("rhoaR: Frame:", i, "/", len(rhoa))
    ax.clear()
    Ert.show(ertData, values=rhoa[i]/rhoa[0], axes=ax, scheme='DipoleDipole', cMin=0.5, cMax=2, 
             colorBar=0, cmap='bwr')
    plt.tight_layout()
createAnimation(fig, animate, int(len(rhoa)), 92, 'rhoaR')   

pg.wait()


