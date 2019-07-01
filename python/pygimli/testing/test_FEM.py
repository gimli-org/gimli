#!/usr/bin/env python
# -*- coding: utf-8 -*-

# write a correct test!
import unittest

import pygimli as pg
import numpy as np

class TestFiniteElementBasics(unittest.TestCase):
    def test_Poisson(self):
        """
        \d² u / d x² + f = 0

        u = sin(x) ## maybe a good test case for harmonic (spectral) base functions TODO
        f = sin(x)

        u = -cos(x)
        f = -cos(x)

        ### the following are already tested in Helmholtz
        u = x 
        f = 0

        u = x*x
        f = -2
        """
        pass 

    def test_Helmholtz(self):
        """
            d² u / d x² + k u + f = 0
            k = 2
            a) P1(exact)
                u = x
                f = -2x
            b) P2(exact)
                u = x*x
                f = -(2 + 2x*x)
        """
        h = np.pi/2 / 21
        x = np.arange(0.0, np.pi/2, h)
        mesh = pg.createGrid(x)

        ### test a)
        k = 2.0
        u = lambda _x: _x
        f = lambda _x: -(k * u(_x))

        x = pg.x(mesh)
        dirichletBC = [[1, u(min(x))], [2, u(max(x))]]
        uFEM = pg.solve(mesh, a=1, b=k, f=f(x), bc={'Dirichlet': dirichletBC})
        np.testing.assert_allclose(uFEM, u(x))

        ### test b)
        u = lambda _x: _x * _x
        f = lambda _x: -(2. + k *u(_x))

        mesh = mesh.createP2()
        x = pg.x(mesh)
        dirichletBC = [[1, u(min(x))], [2, u(max(x))]]
        uFEM = pg.solve(mesh, a=1, b=k, f=f(x), bc={'Dirichlet': dirichletBC})
        np.testing.assert_allclose(uFEM, u(x), atol=1e-6)

        # pg.plt.plot(x, uFEM, '.')
        # pg.plt.plot(pg.sort(x), u(pg.sort(x)))
        # pg.wait()
    
    def test_Neumann(self):
        """
        """
        def _test_(mesh, show=False):
            vTest = 0.1
            u = pg.solve(mesh, a=1, f=0,
                         bc={'Node': [mesh.findNearestNode([0.0, 0.0]), 0.],
                             'Neumann': [[1, -vTest], [2, vTest]]}, verbose=0)

            if show:
                if mesh.dim() == 1:
                    pg.plt.plot(pg.x(mesh), u)
                elif mesh.dim() == 2:
                    pg.show(grid, pg.abs(v))
                    pg.show(grid, v, showMesh=1)
                pg.wait()
            
            v = pg.solver.grad(mesh, u)
            np.testing.assert_allclose(pg.abs(v), np.ones(mesh.cellCount())*vTest)
            return v

        _test_(pg.createGrid(x=np.linspace(-2, 1, 11)), show=False) #1D
        _test_(pg.createGrid(x=np.linspace(-2, 2, 41), y=np.linspace(0, 1, 21))) #2D reg
        _test_(pg.createGrid(x=np.linspace(-0.04, 0.01, 21), y=np.linspace(-0.4, 0, 21))) #2D scaled
        _test_(pg.createGrid(x=np.linspace(-2, 1, 11), y=np.linspace( 0, 1, 11),
                             z=np.linspace( 0, 1, 11))) #3D

        grid = pg.createGrid(x=np.linspace(-2, 2, 41), y=np.linspace(0, 1, 21))
        grid.rotate([0, 0, np.pi/4])
        v = _test_(grid) #2D reg - rotated

        #TODO 2D, Tri, 3D Tet

    def test_Dirichlet(self):
        """
        """
        def _testP1_(mesh, show=False):
            """ Laplace u = 0 solves u = x for u(r=0)=0 and u(r=1)=1
                Test for u == exact x for P1 base functions
            """
            u = pg.solve(mesh, a=1, b=0, f=0, 
                         bc={'Dirichlet': [[1, 0], [2, 1]]})

            if show:
                if mesh.dim()==1:    
                    pg.plt.plot(pg.x(mesh), u)
                    pg.wait()
                elif mesh.dim()==2:
                    pg.show(mesh, u, label='u')
                    pg.wait()

            xMin = mesh.xmin()
            xSpan = (mesh.xmax() - xMin)
            np.testing.assert_allclose(u, (pg.x(mesh)-xMin) / xSpan)
            return u

        def _testP2_(mesh, show=False):
            """ Laplace u = 2 solves u = x² for u(r=0)=0 and u(r=1)=1
                Test for u == exact x² for P2 base functions
            """
            meshp2 = mesh.createP2()
            u = pg.solve(meshp2, f=-2, bc={'Dirichlet': [[1, 0], [2, 1]]})

            # find test pos different from node pos
            meshTests = mesh.createH2()
            meshTests = meshTests.createH2()

            c = [c.center() for c in meshTests.cells()]
            startPos = meshTests.node(0).pos()

            if mesh.dim() == 2:
                c = [b.center() for b in meshTests.boundaries(meshTests.boundaryMarkers()==4)]

            c.sort(key=lambda c_: c_.distance(startPos))
            ui = pg.interpolate(meshp2, u, c)
            xi = pg.utils.cumDist(c) + startPos.distance(c[0])       

            if show:
                pg.plt.plot(xi, ui)
                pg.plt.plot(xi, xi**2)
                pg.wait()
            
            np.testing.assert_allclose(ui, xi**2)
     
        _testP1_(pg.createGrid(x=np.linspace(0, 1, 11)), show=False) #1D
        _testP1_(pg.createGrid(x=np.linspace(-2, 1, 11), 
                               y=np.linspace(0, 1, 11))) #2D reg quad
        _testP1_(pg.createGrid(x=np.linspace(-0.04, 0.01, 11), 
                               y=np.linspace(-0.4, 0, 11))) #2D scaled
        _testP1_(pg.createGrid(x=np.linspace(-2, 1, 11), 
                               y=np.linspace( 0, 1, 11),
                               z=np.linspace( 0, 1, 11))) #3D

        mesh = pg.meshtools.createMesh(pg.meshtools.createWorld(start=[4, -4], end=[6, -6], worldMarker=0), area=0.1)
        mesh.setBoundaryMarkers(np.array([0,1,3,2,4])[mesh.boundaryMarkers()])
        _testP1_(mesh, show=False) #2D tri

        grid = pg.createGrid(x=np.linspace(-2, 2, 11), y=np.linspace(0, 1, 11))
        grid.rotate([0, 0, np.pi/4])
        #can't find proper test for this rotated case
        #v = _testP1_(grid, show=False) #2D reg - rotated
        
        # P2 tests
        mesh = pg.meshtools.createMesh(pg.meshtools.createWorld(start=[0, -4], end=[1, -6], worldMarker=0), area=0.1)
        mesh.setBoundaryMarkers(np.array([0,1,3,2,4])[mesh.boundaryMarkers()])
        _testP2_(mesh, show=False) #2D tri
        _testP2_(pg.createGrid(x=np.linspace(0, 1, 11)), show=0) #1D
        _testP2_(pg.createGrid(x=np.linspace(0, 1, 11), y=np.linspace(0, 1, 11)), show=0) #2D reg quad
        grid = pg.createGrid(x=np.linspace(0, 1, 11), y=np.linspace(0, 1, 11))
        grid.rotate([0, 0, np.pi/4])
        v = _testP2_(grid, show=False) #2D reg - rotated


        #TODO 3D Tet

if __name__ == '__main__':
    
    # test = TestFiniteElementBasics()
    # test.test_Neumann()
    # test.test_Dirichlet()

    unittest.main()