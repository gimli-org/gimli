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
    block = pt.createRectangle(start=[-6, -3.5], end=[6, -6.0], marker=4)
    
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

savefigs_ = 1
savefig(mesh, plc, mesh.cellMarker(), 'Marker', 'marker')
savefig(mesh, plc, mesh.cellMarker(), 'Marker', 'mesh', 1)
savefig(mesh, plc, poro, 'Porosity', 'poro')
savefig(mesh, plc, perm, 'Permeabilty [m$^2$]', 'perm')
savefig(mesh, plc, hydr, 'Hydraulic conductivity [m$/$s]', 'hydr')
savefig(mesh, plc, vP, 'Seismic velocity [m$/$s]', 'vP')
savefig(mesh, plc, dens, 'Density [kg$/$m$^3$]', 'dens')

grav = 1
if grav:
    Grav, densBlock = gravimetry.calcInvBlock(mesh, dens, 'gravInv.pdf')
    savefig(mesh, plc, densBlock, 'Delta Density [kg$/$m$^3$]', 'ddens')
    savefig(Grav.fop.regionManager().paraDomain(), plc,
            Grav.inv.model(), 'Delta Density [kg$/$m$^3$]', 'densGrav')

seis = 0 # geht noch nicht
if seis:
    seismics.calcSeismics(mesh, vP)

# vel brauchen wir ab hier immer
visc = 1./perm
vel = fluidFlow.calcStokesVelocity(mesh, visc, velBoundary, preBoundary)
print("vel:", swatch.duration(True))

showVel = 1
if showVel:
    axVel, _ = pg.show(mesh, data=pg.cellDataToPointData(mesh,
                              np.sqrt(vel[:,0]*vel[:,0] + vel[:,1]*vel[:,1])),
                       logScale=0, colorBar=1,
                       label='|Velocity| in m/s', hold=1)

    pg.viewer.showBoundaryNorm(mesh, velBoundary, color='red', axes=axVel)
    pg.show(plc, axes=axVel)
    
    meshC, _ = createTestWorld1(maxArea=1, verbose=0)
    pg.show(mesh, data=vel, axes=axVel, coarseMesh=meshC, savefig='velocity.pdf')
    

################## TIMELAPSE Starts here 
timeSteps = np.linspace(0, 50, 20) 


################## Concentration ######################
conc = fluidFlow.calcConcentration(mesh, vel, timeSteps, 
                                   injectPos=[-18., -6], scale=1, peclet=50)
print("conc:", swatch.duration(True))


################## TIMELAPSE ERT           ######################
ertMesh, _ = createTestWorld1(maxArea=0.1, verbose=0)

resis, rhoa, ert, ertData = ert.calcApparentResistivities(mesh, ertMesh, poro,
                                          rhoBrine=1./(1./20. + abs(1.*conc)))

# ert inversion als nächstes

pg.wait()


































print("rhoa:", swatch.duration(True))

ertModels, meshERT = calcERT(ert, ertScheme, rhoa)
print("ert:", swatch.duration(True))

vP = velocityVp(porosity, mesh=mesh)
print("vp:", swatch.duration(True))




dpi=92
fig = None

orientation = 'horizontal'

if vis:

    if mp4:
        fig = plt.figure(facecolor='white', figsize=(2*800/dpi, 2*490/dpi), dpi=dpi)  
    else:
        fig = plt.figure() 
    
    axPor = fig.add_subplot(4,4,1)
    axPer = fig.add_subplot(4,4,2)
    axDen = fig.add_subplot(4,4,3)
    axVp  = fig.add_subplot(4,4,4)
    
    axVel = fig.add_subplot(4,4,5)
    axCon = fig.add_subplot(4,4,6)
    axDDe = fig.add_subplot(4,4,7)
    axRes = fig.add_subplot(4,4,8)
    
    axGra = fig.add_subplot(4,4,11)
    axReA = fig.add_subplot(4,4,12)
    
    axERT = fig.add_subplot(4,4,15)
    axERR = fig.add_subplot(4,4,16)
        
    # ** Porosity **
    axPor, _= show(mesh, data=porosity, colorBar=1, 
                orientation=orientation, label='Porosity', axes=axPor)
    show(mesh, axes=axPor)
    #axPor.figure.savefig("poro.pdf",bbox_inches='tight')
    
    # ** Permeabilty **
    axPer,_ = show(mesh, data=permeabilty, colorBar=1, 
         orientation=orientation, label='Permeabilty [m$^2$]', axes=axPer)
    show(mesh, axes=axPer)
    #axPer.figure.savefig("perm.pdf",bbox_inches='tight')
    
    
    axPer,_ = show(mesh, data=hydCond, colorBar=1, 
         orientation=orientation, label='Hydraulic conductivity [m/s]', axes=axPer)
    show(mesh, axes=axPer)
    #axPer.figure.savefig("hydCond.pdf",bbox_inches='tight')
    
    # ** Density **
    show(mesh, data=dens0, colorBar=1, 
         orientation=orientation, label='Density [kg/m$^3$]', axes=axDen)
    show(mesh, axes=axPer)
    
    # ** Velocity abs **
    axVel, cbar = show(mesh, data=pg.cellDataToPointData(mesh,
                                np.sqrt(vel[:,0]*vel[:,0] + vel[:,1]*vel[:,1])),
                       logScale=0, colorBar=1,
                       orientation=orientation, label='|Velocity| in m/s',
                       axes=axVel
                       )

    # ** Velocity vector **
    pg.mplviewer.drawMeshBoundaries(axVel, mesh, fitView=True, hideMesh=1)
    pg.viewer.showBoundaryNorm(mesh, velBoundary, color='red', axes=axVel)
    
    meshC, tmp, tmp, tmp=modelFkt(maxArea=1, verbose=0)
    show(mesh, data=vel, axes=axVel, coarseMesh=meshC)
    
    #axVel.figure.savefig("velocity.pdf",bbox_inches='tight')
    
    # ** vP **
    show(mesh, data=vP, colorBar=1, 
         orientation=orientation, label='Vp m/s', axes=axVp)
    show(mesh, axes=axVp)
    
    # Prepare time lapse figures
    # ** Concentration **
    gciCon= pg.mplviewer.drawModel(axCon, mesh, data=conc[1],
                                   cMin=0, cMax=0.1, 
                                   logScale=False
                                   #cmap='b2r'
                                   )
    cbar = createColorbar(gciCon, orientation=orientation, label='Concentration')
    
    gciDDe = pg.mplviewer.drawModel(axDDe, mesh, data=dDens[1],
                                    cMin=0, cMax=20, 
                                    #cmap='b2r'
                                   )
    cbar = createColorbar(gciDDe, orientation=orientation, label='Delta density in kg/m$^3$')
    
    # ** Resistivity (model) **
    gciRes = pg.mplviewer.drawModel(axRes, meshERT_FOP, 
                                    data=resistivities[0],
                                    cMin=20, cMax=700,
                                    )
    cbar = createColorbar(gciRes, orientation=orientation, label='Resistivity')
    axRes.set_xlim((-20, 20))
    axRes.set_ylim((-14, 00))
    
    # ** Apparent resistivities (data) **
    gciARes = ert.show(ertScheme, values=rhoa[0], axes=axReA, 
                       scheme='DipoleDipole',
                       #cMin=100, cMax=300, 
                       orientation=orientation)

    # ** ERT (model) **
    gciERT = pg.mplviewer.drawModel(axERT, meshERT, 
                                    data=ertModels[0],
                                    cMin=20, cMax=700)
    cbar = createColorbar(gciERT, orientation=orientation, label='Resistivity')
    # ** ERT ratio (model) **
    gciERR = pg.mplviewer.drawModel(axERR, meshERT, 
                                    data=ertModels[0]/ertModels[0],
                                    cMin=1/4, cMax=4, cmap='b2r')
    cbar = createColorbar(gciERR, orientation=orientation, label='Ratio')


def animate(i):
    tic = time.time()
        
    axGra.clear()
    axGra.plot(pg.x(gravPoints), dz[i])
    axGra.plot(pg.x(gravPoints), pg.y(gravPoints), 'v', color='black')
    axGra.set_ylabel('Grav in mGal')
    axGra.set_xlim((-20, 20))
    axGra.set_ylim((0, 0.001))
    axGra.grid()
    
    axReA.clear()
    ert.show(ertScheme, values=rhoa[i], axes=axReA, scheme='DipoleDipole',
             #cMin=100, cMax=300, 
             colorBar=0)
     
    if vis:
        
        pg.mplviewer.setMappableData(gciCon, conc[i], 
                                     #cMin=0, cMax=0.03,
                                     logScale=False)
        gciCon.set_clim(0, 0.1)
        pg.mplviewer.setMappableData(gciDDe, dDens[i], 
                                     cMin=0, cMax=20,
                                     logScale=False)
        pg.mplviewer.setMappableData(gciRes, 
                                     resistivities[i],
                                     cMin=20, cMax=700,
                                     logScale=True)
        pg.mplviewer.setMappableData(gciERT, 
                                     ertModels[i+1],
                                     cMin=20, cMax=700,
                                     logScale=True)
        pg.mplviewer.setMappableData(gciERR, 
                                     ertModels[i+1]/ertModels[0],
                                     cMin=1/4, cMax=4,
                                     logScale=True)
        
    print(i, round(time.time()-tic, 2),
          "t=", round(i *(times[1]-times[0]),1),
          "dt:", round(times[1]-times[0],1),
          "sum mass:", round(sum(dDens[i]*mesh.cellSizes()),1), "kg/m "
          "sum:", sum(conc[i]),
          "dsum:", (sum(conc[i])-sum(conc[i-1])),
          )
    if mp4 or vis:
        pass
        plt.pause(0.001)

#animate(50)
#plt.show()
#for i in range(1, len(times)*2-1):
    #animate(i)

anim = animation.FuncAnimation(fig, animate,
                               frames=int(len(conc)),
                               interval=1)#, blit=True)

solutionName = createCacheName('all', mesh, times)+ "-" + str(ertScheme.size())

if mp4:
    anim.save(solutionName + ".mp4", writer=None, fps=20, dpi=dpi, codec=None,
          bitrate=24*1024, extra_args=None, metadata=None,
          extra_anim=None, savefig_kwargs=None)

pg.showNow()
