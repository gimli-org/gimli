#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys

import numpy as np
import matplotlib.pyplot as plt

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
    fig.tight_layout()

def testColorbar():

    grid = pg.createGrid(x=np.linspace(10., 110., 11)-5, y=np.linspace(0., 20, 2))

    fig, axs = plt.subplots(nrows=3, ncols=3, figsize=((10,6)))
    ax, cbar = pg.show(grid, data=pg.x(grid.cellCenter()), label='log x',
                       ax=axs[0][0], showMesh=True, cMap='Paired', logScale=True)
    ax, cbar = pg.show(grid, data=pg.x(grid.cellCenter())+5., label='log x',
                       ax=axs[1][0], showMesh=True, cMap='Paired', logScale=True)
    pg.viewer.mpl.setMappableData(cbar.mappable, pg.x(grid.cellCenter()))

    ax, cbar = pg.show(grid, data=pg.x(grid.cellCenter()), logScale=True,
                       ax=axs[2][0], showMesh=True, cMap='Paired')
    ###########################################################################
    ax, cbar = pg.show(grid, data=pg.x(grid.cellCenter()), label='lin x',
                       ax=axs[0][1], logScale=False, cMap='Paired')
    ax, cbar = pg.show(grid, data=pg.x(grid.cellCenter())+5, label='lin x',
                       ax=axs[1][1], logScale=False, cMap='Paired')
    pg.viewer.mpl.setMappableData(cbar.mappable, pg.x(grid.cellCenter()))
    ###########################################################################
    grid = pg.createGrid(x=np.linspace(10., 110., 11)-25, y=np.linspace(0., 20, 2))
    ax, cbar = pg.show(grid, data=pg.x(grid.cellCenter()), label='log with neg. x',
                       ax=axs[0][2], showMesh=True, cMap='Paired', logScale=True)
    ax, cbar = pg.show(grid, data=pg.x(grid.cellCenter())-45, label='log with neg. x',
                       ax=axs[1][2], showMesh=True, cMap='Paired', logScale=True)
    pg.viewer.mpl.setMappableData(cbar.mappable, pg.x(grid.cellCenter()))
    ax.figure.tight_layout()

    pg.show(grid, pg.x(grid.cellCenter()), tri=True, shading='gouraud',
            cMap='Spectral_r', logScale=False, cMin=0.01, cMax=10,
            levels=[10, 55, 100],
            orientation="vertical",
            colorBar=True)


    pg.show(grid, pg.x(grid.cellCenter()), tri=True, shading='gouraud',
            cMap='Spectral_r', logScale=False, cMin=0.01, cMax=10,
            levels=[10, 55, 100],
            orientation="horizontal",
            colorBar=True)


def testColRange():
    n = 5
    mesh = pg.createGrid(n, n)

    fig, ax = plt.subplots(2,3, figsize=(8,8))

    ax[0][0].set_title("pyGIMLi (CellData)")
    data = np.arange(mesh.cellCount())
    pg.show(mesh, data, ax=ax[0][0], label='default (nLevs=5, nCols=256)')
    pg.show(mesh, data, ax=ax[1][0], nCols=4, nLevs=5, label="nLevs=5 nCols=4, \n(aka. poor man's contour)")

    ax[0][1].set_title("pyGIMLi (CellData) shading='gouraud'")
    data = np.arange(mesh.cellCount())
    pg.show(mesh, data, ax=ax[0][1], shading='gouraud', label='default (nLevs=5, nCols=256)')
    pg.show(mesh, data, ax=ax[1][1], shading='gouraud', nCols=4, nLevs=5, label="nLevs=5 nCols=4, \n(aka. poor man's contour)")

    ax[0][2].set_title("pyGIMLi (NodeData)")
    data = np.arange(mesh.nodeCount())
    pg.show(mesh, data, ax=ax[0][2],
            label='default (nLevs=5, nCols=4)')
    pg.show(mesh, data, ax=ax[1][2], nCols=12, nLevs=4,
            label='nLevs=4, nCols=12')

    for a in ax.flatten():
        a.yaxis.set_visible(False)
        a.xaxis.set_visible(False)

    fig.tight_layout()

def testCBarLevels():
    """
    Expectations
    ------------
    axs[0, 0]: show regions with plc
        Show needs to deliver the regions with Set3 colormap. Each tick on the
        colobar should be in the middle of the related color section.

    axs[1, 0]: show regions with mesh
        really the same thing as on axs[0, 0] but with mesh.

    ax[0, 1]: show mesh with cell data
        if nLevs is given i would expect that the colormap then is levelled.
        currently that is not the fact. but at least its the full range. labels
        need to be at begin/end of each color section.

    ax[1, 1]: show mesh with node data
        the colorbar range misses parts of its full range. labels need to be
        at begin/end of each color section.
    """
    # create a geometry
    world = mt.createWorld(
        start=[-10, 0], end=[10, -16], layers=[-8], worldMarker=False
        )

    block = mt.createRectangle(
        start=[-6, -3.5], end=[6, -6.0], marker=4, boundaryMarker=10, area=0.1
        )

    circ = mt.createCircle(pos=[0, -11], marker=5, radius=2, boundaryMarker=11)

    poly = mt.mergePLC([world, block, circ])
    mesh = mt.createMesh(poly, quality=34)

    # create random data
    rhomap = [[1, 10], [2, 4], [4, 20], [5, 8]]
    # map data to cell/node count
    cell_data = pg.solver.parseArgToArray(rhomap, mesh.cellCount(), mesh)
    node_data = mt.cellDataToNodeData(mesh, cell_data)

    # plot everything
    fig, axs = pg.plt.subplots(2, 2, figsize=(20, 10))

    pg.show(poly, ax=axs[0, 0])

    pg.show(mesh, ax=axs[1, 0], markers=True)

    pg.show(mesh, cell_data, ax=axs[0, 1], colorBar=True, nLevs=7)

    pg.show(mesh, node_data, ax=axs[1, 1], colorBar=True, nLevs=7)


def testShowPV():
    m1 = mt.createCube()
    # pg.rc['view3D'] = 'fallback'
    # pg.show(m1)

    pg.rc['view3D'] = 'pyvista'
    #pg.show(m1, notebook=True)

    pg.show(m1, gui=False)

def testCoverage():
    grid = pg.createGrid(10,10)
    cov = pg.y(grid.cellCenters()).array()
    cov[-9:] = 0 # remove first row
    data = pg.Vector(grid.cellCount(), 1.0)
    pg.show(grid, data, coverage=cov)
    
if __name__ == '__main__':
    if len(sys.argv) > 0:
        locals()[sys.argv[1]]()
    else:
        testShowVariants()
        testColorbar()
        testCBarLevels()
        testColRange()
        testCoverage()
    
