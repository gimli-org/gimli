#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys

import unittest
import numpy as np
import pygimli as pg


class TestMeshGenerator(unittest.TestCase):

    def test_meshAccess(self):

        x = [0.5, 1, 2, 3, 42]
        y = [0.5, 1, 2, 3]
        z = [0.5, 1, 3]

        mesh = pg.createGrid(x, y, z)

    def test_meshGenRValueProblem(self):
        """ """
        dx = 100
        x = np.arange(-1000, 1001, dx)
        z = np.arange(-800, 1, dx)
        grid = pg.meshtools.createGrid(x=x, y=x, z=z)
        print(grid)
        print(grid.bb())

    def test_triangle(self):
        plc = pg.meshtools.createRectangle()
        mesh = pg.meshtools.createMesh(plc)
        self.assertEqual(mesh.nodeCount(), 4)
        self.assertEqual(mesh.cellCount(), 2)
        self.assertEqual(mesh.boundaryCount(), 5)

    def test_triangle_MAC(self):
        """There seems to be an issue on Mac for higher quality, which produce
        different meshes on mac vs. Linux/Windows. Needs observation and/or
        clarfication."""
        import pygimli as pg
        import pygimli.meshtools as mt

        plc = mt.createCircle(nSegments=24)
        l = mt.createLine(start=[0, -1], end=[0, -0.1], boundaryMarker=2)
        mesh = mt.createMesh([plc, l], area=0.1, quality=30)
        print(mesh)
        # On Linux and Mac: Mesh: Nodes: 43 Cells: 60 Boundaries: 102
        self.assertEqual(mesh.nodeCount(), 43)

        mesh = mt.createMesh([plc, l], area=0.1, quality=32)
        print(mesh)
        if sys.platform == "darwin":
            # On Mac: Mesh: Nodes: 46 Cells: 66 Boundaries: 111
            self.assertEqual(mesh.nodeCount(), 46)
        else:
            # On Linux Mesh: Nodes: 43 Cells: 60 Boundaries: 102 (same as for quality=30)
            self.assertEqual(mesh.nodeCount(), 43)

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

        mesh = pg.meshtools.createMesh2D(np.linspace(0, 1, 6), np.linspace(0, 1, 3))
        self.assertEqual(mesh.cellCount(), 10.0)

    def test_createMesh3D(self):

        mesh = pg.meshtools.createMesh3D(xDim=5, yDim=3, zDim=2)
        self.assertEqual(mesh.cellCount(), 30.0)

    def test_createPartMesh(self):
        mesh = pg.meshtools.createMesh1D(np.linspace(0, 1, 10))
        self.assertEqual(mesh.cellCount(), 9)

        mesh2 = mesh.createMeshByCellIdx(pg.find(pg.x(mesh.cellCenters()) < 0.5))
        self.assertEqual(mesh2.cellCount(), 4)
        self.assertEqual(mesh2.cellCenters()[-1][0] < 0.5, True)

    def test_MeshCreatePolyList(self):
        pos = [[0, 0], [1, 0], [1, -1], [0, -1]]
        poly = pg.meshtools.createPolygon(pos, isClosed=0)
        mesh = pg.meshtools.createMesh(poly, quality=20, area=0.001)
        self.assertEqual(mesh.nodeCount(), 4)
        self.assertEqual(mesh.cellCount(), 0)
        poly = pg.meshtools.createPolygon(pos, isClosed=1)
        mesh = pg.meshtools.createMesh(poly, quality=0, area=0.0)
        self.assertEqual(mesh.nodeCount(), 4)
        self.assertEqual(mesh.cellCount(), 2)
        self.assertEqual(mesh.boundaryCount(), 5)

    def test_MeshCreateSecNodes(self):
        x = [0, 1, 2, 3, 42]
        y = [0, 1, 2, 3]
        z = [0, 1, 3]

        mesh = pg.createGrid(x, y, z)
        mesh.createSecondaryNodes(n=1)
        self.assertEqual(
            mesh.secondaryNodeCount(),
            mesh.boundaryCount()
            + (len(x) - 1) * len(y) * len(z)
            + (len(y) - 1) * len(z) * len(x)
            + (len(z) - 1) * len(x) * len(y),
        )

    def test_MeshStr(self):
        mesh = pg.createGrid(2, 2, 2)
        print(mesh.node(0))

    def test_MeshDataAccess(self):
        mesh = pg.Mesh()
        a = pg.Vector(10, 1.0)
        b = [pg.Vector(10, 1.0)] * 3
        c = np.array(b).T

        mesh["a"] = a
        mesh["b"] = b
        mesh["v"] = c
        mesh["vs"] = [c, c, c]

        # pg.core.setDeepDebug(True)
        # pg.core.setDeepDebug(False)

        np.testing.assert_array_equal(mesh["a"], a)
        np.testing.assert_array_equal(mesh["b"], b)
        np.testing.assert_array_equal(mesh["v"], c)
        np.testing.assert_array_equal(mesh["vs"], [c, c, c])

        # mesh['c'] = pg.PosList(10, [1.0, 0., 0.0])

    def test_meshBMS(self):
        # text bms version v3 which stores geometry flag
        mesh = pg.Mesh(2, isGeometry=True)

        import tempfile as tmp

        _, fn = tmp.mkstemp()

        mesh.save(fn)
        mesh2 = pg.load(fn + ".bms", verbose=True)

        self.assertEqual(mesh.isGeometry(), mesh2.isGeometry())

    def test_VTK_DataRead(self):
        grid = pg.createGrid(np.arange(4), np.arange(3), np.arange(2))
        cM = np.arange(grid.cellCount())
        grid.setCellMarkers(cM)

        import tempfile as tmp

        _, fn = tmp.mkstemp(suffix=".vtk")

        grid.exportVTK(fn)
        mesh = pg.load(fn)
        np.testing.assert_array_equal(mesh.cellMarkers(), cM)
        np.testing.assert_array_equal(mesh["Marker"], cM)

        mesh = pg.meshtools.readMeshIO(fn)
        np.testing.assert_array_equal(mesh["Marker"], cM)

        fn = pg.getExampleFile("meshes/test_tetgen_dataCol.vtk")
        mesh = pg.load(fn)
        np.testing.assert_array_equal(mesh.cellMarkers(), cM)
        np.testing.assert_array_equal(mesh["Marker"], cM)

        mesh = pg.meshtools.readMeshIO(fn)
        np.testing.assert_array_equal(mesh["Marker"], cM)

        # pg._g('pg import vtk')
        # print(mesh)
        # print(mesh["Marker"])
        # print(mesh.cellMarkers())

        # pg._g('pg import tetgen vtk')
        # mesh = pg.load("grid1.vtk")
        # print(mesh)
        # print(mesh["Marker"])
        # print(mesh.cellMarkers())

        # pg._g('meshio import vtk')
        # print(mesh)
        # print(mesh["Marker"])
        # #print(mesh.cellMarkers())

        # pg._g('meshio import tetgen vtk')
        # mesh = pg.meshtools.readMeshIO("grid1.vtk")
        # print(mesh)
        # print(mesh["Marker"])

    def test_VTK_ExportVTU(self):
        """Test to fix export bug"""
        mesh = pg.createGrid(4, 4, 4)
        mesh.exportBoundaryVTU("bounds.vtu")

    def test_SimpleMeshExport(self):

        mesh = pg.createGrid(3, 3)
        verts = mesh.positions()
        cellIds = [c.ids() for c in mesh.cells()]

        mesh2 = pg.Mesh(2)
        [mesh2.createNode(v) for v in verts]
        [mesh2.createCell(c) for c in cellIds]

        np.testing.assert_array_equal(mesh2.nodeCount(), mesh.nodeCount())
        np.testing.assert_array_equal(mesh2.cellCount(), mesh.cellCount())


if __name__ == "__main__":
    # pg.setDeepDebug(1)

    # t = TestMeshGenerator()
    # t.test_MeshDataAccess()
    # sys.exit()
    # t.test_meshAccess()
    # t.test_createGrid()
    # exit()

    unittest.main()
