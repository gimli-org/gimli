#!/usr/bin/env python
# -*- coding: utf-8 -*-
import unittest
import numpy as np
import pygimli as pg


class TestMeshGenerator(unittest.TestCase):

    def test_meshAccess(self):

        x = [.5, 1, 2, 3, 42]
        y = [.5, 1, 2, 3]
        z = [.5, 1, 3]

        mesh = pg.createGrid(x, y, z)

    def test_triangle(self):
        plc = pg.meshtools.createRectangle()
        mesh = pg.meshtools.createMesh(plc)
        self.assertEqual(mesh.nodeCount(), 4)
        self.assertEqual(mesh.cellCount(), 2)
        self.assertEqual(mesh.boundaryCount(), 5)

    def test_createGrid(self):
        mesh = pg.createGrid(3)
        self.assertEqual(mesh.xmax(), 2.0)
        mesh = pg.createGrid(3, 3)
        self.assertEqual(mesh.ymax(), 2.0)
        mesh = pg.createGrid(3, 3, 3)
        self.assertEqual(mesh.zmax(), 2.0)
        # mesh = pg.meshtools.createMesh1D(10, 1)
        # print(mesh)
        # self.assertEqual(mesh.cellCount(), 10.0)

    def test_createMesh1D(self):

        mesh = pg.meshtools.createMesh1D(10, 1)
        self.assertEqual(mesh.cellCount(), 10.0)
        self.assertEqual(mesh.xmax(), 10.0)

        mesh = pg.meshtools.createMesh1D(nCells=10)
        self.assertEqual(mesh.cellCount(), 10.0)

        mesh = pg.meshtools.createMesh1D(nCells=5, nProperties=2)
        self.assertEqual(mesh.cellCount(), 10.0)

        mesh = pg.meshtools.createMesh1D(5, 2)
        self.assertEqual(mesh.cellCount(), 10.0)

        mesh = pg.meshtools.createMesh1D(10)
        self.assertEqual(mesh.cellCount(), 10.0)

    def test_createMesh1DBlock(self):

        mesh = pg.meshtools.createMesh1DBlock(nLayers=5)
        self.assertEqual(mesh.cellCount(), 9.0)

        mesh = pg.meshtools.createMesh1DBlock(5)
        self.assertEqual(mesh.cellCount(), 9.0)

        mesh = pg.meshtools.createMesh1DBlock(5, 1)
        self.assertEqual(mesh.cellCount(), 9.0)

        mesh = pg.meshtools.createMesh1DBlock(nLayers=4, nProperties=2)
        self.assertEqual(mesh.cellCount(), 11.0)

        mesh = pg.meshtools.createMesh1DBlock(4, 2)
        self.assertEqual(mesh.cellCount(), 11.0)

    def test_createMesh2D(self):

        mesh = pg.meshtools.createMesh2D(xDim=5, yDim=2)
        self.assertEqual(mesh.cellCount(), 10.0)

        mesh = pg.meshtools.createMesh2D(5, 2)
        self.assertEqual(mesh.cellCount(), 10.0)

        mesh = pg.meshtools.createMesh2D(np.linspace(0, 1, 6),np.linspace(0, 1, 3))
        self.assertEqual(mesh.cellCount(), 10.0)

    def test_createMesh3D(self):

        mesh = pg.meshtools.createMesh3D(xDim=5, yDim=3, zDim=2)
        self.assertEqual(mesh.cellCount(), 30.0)

    def test_createPartMesh(self):
        mesh = pg.meshtools.createMesh1D(np.linspace(0, 1, 10))
        self.assertEqual(mesh.cellCount(), 9)

        mesh2 = mesh.createMeshByCellIdx(
            pg.find(pg.x(mesh.cellCenters()) < 0.5))
        self.assertEqual(mesh2.cellCount(), 4)
        self.assertEqual(mesh2.cellCenters()[-1][0] < 0.5, True)

    def test_MeshCreatePolyList(self):
        pos = [[0, 0], [1, 0], [1, -1], [0, -1]]
        poly = pg.meshtools.createPolygon(pos, isClosed=0)
        mesh = pg.meshtools.createMesh(poly, quality=20, area=0.001)
        self.assertEqual(mesh.nodeCount(), 4)
        self.assertEqual(mesh.cellCount(), 0)
        poly = pg.meshtools.createPolygon(pos, isClosed=1)
        mesh = pg.meshtools.createMesh(poly, quality=0, area=0.)
        self.assertEqual(mesh.nodeCount(), 4)
        self.assertEqual(mesh.cellCount(), 2)
        self.assertEqual(mesh.boundaryCount(), 5)

    def test_MeshCreateSecNodes(self):
        x = [0, 1, 2, 3, 42]
        y = [0, 1, 2, 3]
        z = [0, 1, 3]

        mesh = pg.createGrid(x, y, z)
        mesh.createSecondaryNodes(n=1)
        self.assertEqual(mesh.secondaryNodeCount(), mesh.boundaryCount() + \
                                                    (len(x)-1)*len(y)*len(z) + \
                                                    (len(y)-1)*len(z)*len(x) + \
                                                    (len(z)-1)*len(x)*len(y))

    def test_MeshStr(self):
        mesh= pg.createGrid(2,2,2)
        print(mesh.node(0))


if __name__ == '__main__':
    # pg.setDeepDebug(1)
    
    # t = TestMeshGenerator()
    # t.test_meshAccess()
    # t.test_createGrid()
    # exit()

    unittest.main()
