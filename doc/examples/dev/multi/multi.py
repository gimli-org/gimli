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

def savefig(mesh, plc, data=None, label='', out=None, showMesh=0):

    ax = None
    if data is not None:
        ax, _ = pg.show(mesh, data, hold=1, colorBar=1, pad=0.55, label=label)
    if showMesh:
        ax, _ = pg.show(mesh, axes=ax, hold=1)
    if out:
        ax, _ = pg.show(plc, axes=ax)
        
        adjustAxes(ax)
        plt.pause(0.01)
        ax.figure.savefig(out + '.pdf', bbox_inches='tight')

        try:
            print("trying pdf2pdfS ... ")
            os.system('pdf2pdfBB ' + out + '.pdf')
            os.system('pdf2pdfS ' + out + '.pdf')
        except:
            pass
    else:
        ax, _ = pg.show(plc, axes=ax)
    return ax

def adjustAxes(ax):
    ax.set_ylabel('Depth [m]')
    ax.set_xlabel('$x$ [m]')
        
    ticks = ax.yaxis.get_majorticklocs()
    tickLabels = []
    for t in ticks:
        tickLabels.append(str(int(abs(t))))

    ax.set_yticklabels(tickLabels)

def createAnimation(fig, animate, nFrames, dpi, out):
    anim = animation.FuncAnimation(fig, animate,
                                   frames=nFrames,
                                   interval=0.001, repeat=False)
    anim.save(out + ".mp4", writer=None, fps=20, dpi=dpi, codec=None,
              bitrate=24*1024, extra_args=None, metadata=None,
              extra_anim=None, savefig_kwargs=None)
    try:
        print("Create frames ... ")
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
    
    cbar = pg.mplviewer.createColorbar(gci, label=label, pad=0.55)
    ax.set_ylabel('Depth [m]')
    ax.set_xlabel('$x$ [m]')
        
    ticks = ax.yaxis.get_majorticklocs()
    tickLabels = []
    for t in ticks:
        tickLabels.append(str(int(abs(t))))

    ax.set_yticklabels(tickLabels)
    
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
                               boundaryMarker=10, area=0.1)
    
    plc = pt.mergePLC([layer1, layer2, layer3, block])
    
    mesh = createMesh(plc, quality=33, area=maxArea,
                      smooth=[1,10], verbose=verbose)
   
    return mesh, plc
    
swatch = pg.core.Stopwatch(True)

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

vP   = seismics.velocityVp(poro, mesh=mesh)

seis = 1
if seis:
    seismics.calcSeismics(mesh, vP)

pg.wait()


