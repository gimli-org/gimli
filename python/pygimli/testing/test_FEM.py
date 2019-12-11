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
        def _test_(mesh, show=False):
            """
                \Laplace u = 0
                du/dn = -0.1 (xMin)
                du/dn =  0.1 (xMax)
            """
            vTest = 0.1
            u = pg.solve(mesh, a=1, f=0,
                         bc={'Node': [mesh.findNearestNode([0.0, 0.0]), 0.],
                             'Neumann': [[1, -vTest], [2, vTest]]}, verbose=0)
            v = pg.solver.grad(mesh, u)

            if show:
                if mesh.dim() == 1:
                    pg.plt.plot(pg.x(mesh), u)
                elif mesh.dim() == 2:
                    pg.show(mesh, u, label='u')
                    ax, _ = pg.show(mesh, pg.abs(v), label='abs(grad(u))')
                    pg.show(mesh, v, showMesh=1, ax=ax)
                pg.wait()

            # print(pg.x(mesh))
            # print(pg.x(mesh)*vTest)
            # print(u)

            # print(min(u), max(u))
            # print(min(pg.x(mesh)*vTest), max(pg.x(mesh)*vTest))
            # print(min(pg.x(mesh)*vTest-u), max(pg.x(mesh)*vTest-u))

            # for +1.0 .. see https://github.com/numpy/numpy/issues/13801
            np.testing.assert_allclose(1.0 + u, 1.0 + pg.x(mesh)*vTest)

            np.testing.assert_allclose(pg.abs(v),
                                       np.ones(mesh.cellCount())*vTest)
            return v

        # 1D
        _test_(pg.createGrid(x=np.linspace(-1, 1, 11)), show=False)
        # 1D scaled
        _test_(pg.createGrid(x=np.linspace(-2, 2, 11)))

        # 2D grid
        _test_(pg.createGrid(x=np.linspace(-2, 2, 21),
                             y=np.linspace(0, 1, 11)), show=False)
        # 2D scaled
        _test_(pg.createGrid(x=np.linspace(-0.04, 0.04, 21),
                             y=np.linspace(-0.4, 0, 21)))

        # 3D grid
        _test_(pg.createGrid(x=np.linspace(-2, 2, 11), y=np.linspace(0, 1, 11),
                             z=np.linspace( 0, 1, 11)))

        # 2D grid rotated -- need better test
        # grid = pg.createGrid(x=np.linspace(-2, 2, 41), y=np.linspace(-2, 2, 21))
        # grid.rotate([0, 0, np.pi])
        # _test_(grid, show=True)

        # 2D tri
        mesh = pg.meshtools.createMesh(pg.meshtools.createWorld(start=[-1, -4],
                                                                end=[1, -6], worldMarker=0), area=0.1)
        mesh.setBoundaryMarkers(np.array([0,1,3,2,4])[mesh.boundaryMarkers()])
        _test_(mesh, show=False)

        #TODO 2D, Tri, 3D Tet

    def test_Dirichlet(self):
        """
        """
        def _testP2_(mesh, show=False):
            """ Laplace u - 2 = 0
                with u(x) = x² + x/(xMax-xMin) + xMin/(xMax-xMin)
                and u(xMin)=0 and u(xMax)=1
                Test for u_h === u(x) for P2 base functions
            """
            meshP2 = mesh.createP2()
            u = pg.solve(meshP2, f=-2, bc={'Dirichlet': [[1, 2], [2, 2]]},)

            xMin = mesh.xMin()
            xSpan = (mesh.xMax() - xMin)

            if show:
                if mesh.dim()==1:
                    pg.plt.figure()
                    x = pg.x(mesh)
                    ix = np.argsort(x)
                    pg.plt.plot(x[ix], x[ix]**2 - (xSpan/2)**2 + 2 )
                elif mesh.dim() > 1:
                    pg.show(meshP2, u, label='u = x**2')

            np.testing.assert_allclose(u, pg.x(meshP2)**2 - (xSpan/2)**2 +2)


        def _testP1_(mesh, show=False, followP2=True):
            """ Laplace u = 0
                with u = x and u(x=min)=0 and u(x=max)=1
                Test for u == exact x for P1 base functions
            """
            u = pg.solve(mesh, a=1, b=0, f=0,
                         bc={'Dirichlet': [[1, 0], [2, 1]]})

            if show:
                if mesh.dim()==1:
                    pg.plt.figure()
                    x = pg.x(mesh)
                    ix = np.argsort(x)
                    pg.plt.plot(x[ix], u[ix])
                elif mesh.dim() > 1:
                    pg.show(mesh, u, label='u')

            xMin = mesh.xMin()
            xSpan = (mesh.xMax() - xMin)
            np.testing.assert_allclose(u, (pg.x(mesh)-xMin) / xSpan)

            if followP2:
                _testP2_(mesh, show)
            return u


        # 1D
        _testP1_(pg.createGrid(x=np.linspace(-1, 1, 11)), show=False)
        _testP1_(pg.createGrid(x=np.linspace(-3, 3, 11)), show=False)

        # 2D reg quad
        _testP1_(pg.createGrid(x=np.linspace(-2, 2, 11),
                               y=np.linspace(0, 1, 11)), show=False)
        # 2D reg quad scale
        _testP1_(pg.createGrid(x=np.linspace(-0.04, 0.04, 11),
                               y=np.linspace(-0.4, 0, 11)))
        # 3D reg quad
        _testP1_(pg.createGrid(x=np.linspace(-2, 2, 11),
                               y=np.linspace( 0, 1, 11),
                               z=np.linspace( 0, 1, 11)),
                               followP2=False
                               )
                               # check why P2 fails here,

        # 2D tri
        mesh = pg.meshtools.createMesh(pg.meshtools.createWorld(start=[-1, -4],
                                                                end=[1, -6], worldMarker=0), area=0.1)
        mesh.setBoundaryMarkers(np.array([0,1,3,2,4])[mesh.boundaryMarkers()])
        _testP1_(mesh, show=False)

        # 3D prism
        mesh = pg.meshtools.extrudeMesh(mesh, np.linspace(0, 1, 11))
        _testP1_(mesh, show=False)

        # 3D tet
        mesh = pg.meshtools.createMesh(pg.meshtools.createCube(size=[4, 4, 4],
                                                               boundaryMarker=9,
                                                               area=100.1))

        for b in mesh.boundaries(mesh.boundaryMarkers() == 9):

            if b.norm()[0] == -1:
                b.setMarker(1)
            elif b.norm()[0] == +1:
                b.setMarker(2)

        _testP1_(mesh, show=False)
        #pg.wait()

        grid = pg.createGrid(x=np.linspace(-2, 2, 11), y=np.linspace(0, 1, 11))
        grid.rotate([0, 0, np.pi/4])
        #can't find proper test for this rotated case
        #v = _testP1_(grid, show=False) #2D reg - rotated

    def testElementMatrix(self):
        a = pg.core.ElementMatrix()


if __name__ == '__main__':

    # test = TestFiniteElementBasics()
    # test.test_Neumann()
    # test.test_Dirichlet()

    unittest.main()