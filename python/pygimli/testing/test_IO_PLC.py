#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest

import pygimli as pg
import numpy as np


class TestPLCIO(unittest.TestCase):

    def test_io_trinagle(self):
        """
        """
        # create tempfile in most secure manner, only accesible by this process
        # id no execution allowed at all, will be deleted as soon as this
        # process stops
        try:
            import tempfile as tmp
        except ImportError:
            return

        _, name2D = tmp.mkstemp(suffix='.poly')

        # 2D, creating test trinangle poly
        m = pg.Mesh(2)
        nodes = []
        nodes.append(m.createNode([-1, -2]))
        nodes.append(m.createNode([-1, 2]))
        nodes.append(m.createNode([1, 2]))
        nodes.append(m.createNode([1, -2]))
        for i in range(4):
            m.createEdge(nodes[i], nodes[(i + 1) % 4])

        nodes.append(m.createNode([-1, -3]))
        nodes.append(m.createNode([1, -3]))

        m.createEdge(nodes[0], nodes[4])
        m.createEdge(nodes[1], nodes[5])
        m.createEdge(nodes[4], nodes[5])

        m.addRegionMarker([0., 0.], -3, area=-1)
        m.addRegionMarker([0., -2.5], 1, area=42.42)
        pg.meshtools.exportPLC(m, name2D)

        # 2D, test triangle mesh for consistincy
        poly = pg.meshtools.readPLC(name2D)

        np.testing.assert_allclose(poly.regionMarker()[0].array(),
                                   np.array([0., -2.5, 0.]))
        np.testing.assert_allclose(poly.holeMarker()[0].array(),
                                   np.array([0., 0., 0.]))
        np.testing.assert_allclose(np.sort(poly.positions().array()),
                                   np.sort(m.positions().array()))
        np.testing.assert_equal(poly.regionMarker()[0].area(), 42.42)

    def test_io_tetgen(self):
        """
        """
        try:
            import tempfile as tmp
        except ImportError:
            return

        _, name3D = tmp.mkstemp(suffix='.poly')
        print(name3D)

        # 3D, creating test tetgen poly
        minx = -2.0
        maxx = np.pi
        miny = 0.0
        maxy = 23
        minz = -3
        maxz = 1e-4
        cube = pg.Mesh(3)

        ob0 = cube.createNode(minx, miny, minz)
        ob1 = cube.createNode(maxx, miny, minz)
        ob2 = cube.createNode(maxx, maxy, minz)
        ob3 = cube.createNode(minx, maxy, minz)
        ob4 = cube.createNode(minx, miny, maxz)
        ob5 = cube.createNode(maxx, miny, maxz)
        ob6 = cube.createNode(maxx, maxy, maxz)
        ob7 = cube.createNode(minx, maxy, maxz)

        cube.createQuadrangleFace(ob0, ob1, ob2, ob3)
        cube.createQuadrangleFace(ob1, ob2, ob6, ob5)
        cube.createQuadrangleFace(ob0, ob1, ob5, ob4)
        cube.createQuadrangleFace(ob0, ob3, ob7, ob4)
        cube.createQuadrangleFace(ob2, ob3, ob7, ob6)
        cube.createQuadrangleFace(ob4, ob5, ob6, ob7)
        cube.addRegionMarker([0., 0., -2.0], -3, area=-1)
        cube.addRegionMarker([-1.99, 0.001, 1e-6], 1, area=42.42)
        cube.addRegionMarker([-0.99, 1, 1e-7], 1, area=1)
        cube.addRegionMarker([0.99, 11, 1e-8], 1, area=2)
        cube.addRegionMarker([1.99, 22, 1e-9], 1, area=1245535642455)
        pg.meshtools.exportPLC(cube, name3D)

        # 3D, test tetgen plc for consistincy
        poly = pg.meshtools.readPLC(name3D)

        np.testing.assert_allclose(poly.regionMarker()[0].array(),
                                   np.array([-1.99, 0.001, 1e-6]))
        np.testing.assert_allclose(poly.holeMarker()[0].array(),
                                   cube.holeMarker()[0].array())
        np.testing.assert_allclose(np.sort(poly.positions().array()),
                                   np.sort(cube.positions().array()))
        np.testing.assert_equal(poly.regionMarker()[0].area(), 42.42)

if __name__ == '__main__':
    unittest.main()
