#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

#plt.xkcd()

import pygimli as pg
import pygimli.meshtools as mt

def testShowVariants():
    # Create geometry definition for the modelling domain
    world = mt.createWorld(start=[-10, 0], end=[10, -16],
                        layers=[-8], worldMarker=False)

    # Create a heterogeneous block
    block = mt.createRectangle(start=[-6, -3.5], end=[6, -6.0],
                            marker=4,  boundaryMarker=10, area=0.1)

    circ = mt.createCircle(pos=[0, -11], radius=2, boundaryMarker=11, isHole=True)

    # Merge geometrical entities
    geom = world + block + circ
    mesh = mt.createMesh(geom)

    fig, axs = plt.subplots(3,5)

    pg.show(geom, ax=axs[0][0])
    axs[0][0].set_title('plc, (default)')
    pg.show(geom, fillRegion=False, ax=axs[0][1])
    axs[0][1].set_title('plc, fillRegion=False')
    pg.show(geom, showBoundary=False, ax=axs[0][2])
    axs[0][2].set_title('plc, showBoundary=False')
    pg.show(geom, markers=True, ax=axs[0][3])
    axs[0][3].set_title('plc, markers=True')
    pg.show(mesh, ax=axs[0][4], showBoundary=False)
    axs[0][4].set_title('mesh, showBoundary=False')

    pg.show(mesh, ax=axs[1][0])
    axs[1][0].set_title('mesh, (default)')
    pg.show(mesh, mesh.cellMarkers(), label='Cell markers', ax=axs[1][1])
    axs[1][1].set_title('mesh, cells, (default)')
    pg.show(mesh, markers=True, ax=axs[1][2])
    axs[1][2].set_title('mesh, cells, markers=True')
    pg.show(mesh, mesh.cellMarkers(), label='Cell markers', showMesh=True, ax=axs[1][3])
    axs[1][3].set_title('mesh, cells, showMesh=True')
    pg.show(mesh, mesh.cellMarkers(), label='Cell markers', showBoundary=False, ax=axs[1][4])
    axs[1][4].set_title('mesh, cells, showBoundary=False')

    pg.show(mesh, pg.x(mesh), label='Nodes (x)', ax=axs[2][0])
    axs[2][0].set_title('mesh, nodes, (default)')
    pg.show(mesh, pg.x(mesh), label='Nodes (x)', showMesh=True, ax=axs[2][1])
    axs[2][1].set_title('mesh, nodes, showMesh=True')
    pg.show(mesh, pg.x(mesh), label='Nodes (x)', showBoundary=True, ax=axs[2][2])
    axs[2][2].set_title('mesh, nodes, showBoundary=True')
    pg.show(mesh, pg.y(mesh.cellCenters()), label='Cell center (y)',
            tri=True, shading='flat', ax=axs[2][3])
    axs[2][3].set_title('mesh, cells, tri=True, shading=flat')
    pg.show(mesh, pg.y(mesh.cellCenters()), label='Cell center (y)',
            tri=True, shading='gouraud', ax=axs[2][4])
    axs[2][4].set_title('mesh, cells, tri=True, shading=gouraud')
    ##pg.show(mesh, mesh.cellMarker(), label(markers), axs[1][1])
    axs[2][4].figure.tight_layout()

def testColorbar():

    grid = pg.createGrid(x=np.linspace(10., 110., 11)-5, y=np.linspace(0., 20, 2))

    fig, axs = plt.subplots(nrows=3, ncols=3, figsize=((10,6)))
    ax, cbar = pg.show(grid, data=pg.x(grid.cellCenter()), label='log x',
                       ax=axs[0][0], showMesh=True, cMap='Paired', logScale=True)
    ax, cbar = pg.show(grid, data=pg.x(grid.cellCenter())+5., label='log x',
                       ax=axs[1][0], showMesh=True, cMap='Paired', logScale=True)
    pg.mplviewer.setMappableData(cbar.mappable, pg.x(grid.cellCenter()))

    ax, cbar = pg.show(grid, data=pg.x(grid.cellCenter()), logScale=True,
                       ax=axs[2][0], showMesh=True, cMap='Paired')
    ###########################################################################
    ax, cbar = pg.show(grid, data=pg.x(grid.cellCenter()), label='lin x',
                       ax=axs[0][1], logScale=False, cMap='Paired')
    ax, cbar = pg.show(grid, data=pg.x(grid.cellCenter())+5, label='lin x',
                       ax=axs[1][1], logScale=False, cMap='Paired')
    pg.mplviewer.setMappableData(cbar.mappable, pg.x(grid.cellCenter()))
    ###########################################################################
    grid = pg.createGrid(x=np.linspace(10., 110., 11)-25, y=np.linspace(0., 20, 2))
    ax, cbar = pg.show(grid, data=pg.x(grid.cellCenter()), label='log with neg. x',
                       ax=axs[0][2], showMesh=True, cMap='Paired', logScale=True)
    ax, cbar = pg.show(grid, data=pg.x(grid.cellCenter())-45, label='log with neg. x',
                       ax=axs[1][2], showMesh=True, cMap='Paired', logScale=True)
    pg.mplviewer.setMappableData(cbar.mappable, pg.x(grid.cellCenter()))

    ax.figure.tight_layout()


if __name__ == '__main__':
    #testShowVariants()
    testColorbar()
    pg.wait()
