#!/usr/bin/env python
# -*- coding: utf-8 -*-

# write a correct test!
import unittest

import pygimli as pg
import numpy as np

class TestFiniteElementBasics(unittest.TestCase):

    def test_Neumann(self):
        """
        """
        def _test_(mesh):
            vTest = 0.1
            u = pg.solve(mesh, a=1,
                        bc={'Node': [mesh.findNearestNode([0.0, 0.0]), 0.],
                            'Neumann': [[1, -vTest], [2, vTest]]})

            v = pg.solver.grad(mesh, u)
            #print("|v|:", min(pg.abs(v)), max(pg.abs(v)), pg.mean(pg.abs(v)))
            np.testing.assert_allclose(pg.abs(v), np.ones(mesh.cellCount())*vTest)
            return v

        _test_(pg.createGrid(x=np.linspace(-2, 1, 11))) #1D
        _test_(pg.createGrid(x=np.linspace(-2, 2, 41), y=np.linspace(0, 1, 21))) #2D reg
        _test_(pg.createGrid(x=np.linspace(-0.04, 0.01, 21), y=np.linspace(-0.4, 0, 21))) #2D scaled
        _test_(pg.createGrid(x=np.linspace(-2, 1, 11), y=np.linspace( 0, 1, 11),
                             z=np.linspace( 0, 1, 11))) #3D

        grid = pg.createGrid(x=np.linspace(-2, 2, 41), y=np.linspace(0, 1, 21))
        grid.rotate([0, 0, np.pi/4])
        v = _test_(grid) #2D reg - rotated
        #pg.show(grid, pg.abs(v))
        #pg.show(grid, v, showMesh=1)
        #pg.wait()

        #TODO 2D, Tri, 3D Tet

    def test_Dirichlet(self):
        """
        """
        def _test_(mesh, show=False):
            u = pg.solve(mesh, bc={'Dirichlet': [[1, 0], [2, 1]]})

            if show:
                if mesh.dim()==1:    
                    pg.plt.plot(pg.x(mesh), u)
                    pg.wait()
                elif mesh.dim()==2:
                    pg.show(mesh, u, label='u')
                    pg.wait()

            #print("|v|:", min(pg.abs(v)), max(pg.abs(v)), pg.mean(pg.abs(v)))
            xMin = mesh.xmin()
            xSpan = (mesh.xmax() - xMin)
            np.testing.assert_allclose(u, (pg.x(mesh)-xMin)/ xSpan)
            return u

        _test_(pg.createGrid(x=np.linspace(0, 1, 11))) #1D
        _test_(pg.createGrid(x=np.linspace(-2, 1, 41), y=np.linspace(0, 1, 21))) #2D reg quad

        mesh = pg.meshtools.createMesh(pg.meshtools.createWorld(start=[4, -4], end=[6, -6], worldMarker=0), area=0.1)
        mesh.setBoundaryMarkers(np.array([0,1,3,2,4])[mesh.boundaryMarkers()])
        _test_(mesh, show=False) #2D tri

        _test_(pg.createGrid(x=np.linspace(-0.04, 0.01, 21), y=np.linspace(-0.4, 0, 21))) #2D scaled
        _test_(pg.createGrid(x=np.linspace(-2, 1, 11), y=np.linspace( 0, 1, 11),
                             z=np.linspace( 0, 1, 11))) #3D

        # can't find proper test for this rotated case
        # grid = pg.createGrid(x=np.linspace(-2, 2, 41), y=np.linspace(0, 1, 21))
        # grid.rotate([0, 0, np.pi/4])
        # v = _test_(grid, show=True) #2D reg - rotated
        

        #TODO 2D, Tri, 3D Tet

if __name__ == '__main__':
    test = TestFiniteElementBasics()
    test.test_Dirichlet()
    
    #test.test_Neumann()

    #unittest.main()
