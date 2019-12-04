"""Tests for pygimli.meshtools.polytools
"""
import numpy as np
import unittest

import pygimli as pg
import pygimli.meshtools as mt

class TestCreateRectangle(unittest.TestCase):
    def test_region_marker_position_basics(self):
        rect1 = mt.createRectangle(
            start=[0.0, 0.0],
            end=[2.0, -1.0],
            isClosed=True,
            marker=1,
        )
        # by default the region marker should be located at
        # sPos + (ePos - sPos) * 0.2)
        assert rect1.regionMarkers()[0].x() == 0.4
        assert rect1.regionMarkers()[0].y() == -0.2
        assert rect1.regionMarkers()[0].marker() == 1

    def test_region_marker_position_translation_scale(self):
        rect1 = mt.createRectangle(
            pos=[1.0, 2.0],
            size=[2.0, 4],
            isClosed=True,
            marker=20,
        )
        assert rect1.regionMarkers()[0].x() == -0.3
        assert rect1.regionMarkers()[0].y() == 0.3
        assert rect1.regionMarkers()[0].marker() == 20

    def test_region_marker_position_two_ways_v1(self):
        rect1 = mt.createRectangle(
            start=[-2.0, -1.0],
            end=[2.0, 1.0],
            isClosed=True,
            marker=1,
        )

        # scaled version
        rect2 = mt.createRectangle(
            pos=[0.0, 0.0],
            size=[4.0, 2],
            isClosed=True,
            marker=2,
        )
        assert rect1.regionMarkers()[0].x() == -1.2
        assert rect1.regionMarkers()[0].y() == -0.6

        assert rect1.regionMarkers()[0].marker() == 1
        assert rect2.regionMarkers()[0].marker() == 2

        assert rect1.regionMarkers()[0].x() == rect1.regionMarkers()[0].x()
        assert rect2.regionMarkers()[0].x() == rect2.regionMarkers()[0].x()

        assert rect1.regionMarkers()[0].y() == rect1.regionMarkers()[0].y()
        assert rect2.regionMarkers()[0].y() == rect2.regionMarkers()[0].y()

    def test_region_marker_position_two_ways_v2(self):
        rect1 = mt.createRectangle(
            start=[0.0, -1.5],
            end=[4.0, -4.0],
            isClosed=True,
            marker=1,
        )

        # scaled version
        rect2 = mt.createRectangle(
            pos=[2.0, -2.75],
            size=[4.0, 2.5],
            isClosed=True,
            marker=2,
        )
        assert rect1.regionMarkers()[0].x() == 4 * 0.2
        assert rect1.regionMarkers()[0].y() == -1.5 - 2.5 * 0.2

        assert rect1.regionMarkers()[0].marker() == 1
        assert rect2.regionMarkers()[0].marker() == 2

        assert rect1.regionMarkers()[0].x() == rect1.regionMarkers()[0].x()
        assert rect2.regionMarkers()[0].x() == rect2.regionMarkers()[0].x()

        assert rect1.regionMarkers()[0].y() == rect1.regionMarkers()[0].y()
        assert rect2.regionMarkers()[0].y() == rect2.regionMarkers()[0].y()

    def test_region_markerposition_start_end(self):
        rect1 = mt.createRectangle(
            start=[0.0, -1.5],
            end=[4.0, -4.0],
            isClosed=True,
            marker=10,
            markerPosition=[0.0, 0.0],
        )
        assert rect1.regionMarkers()[0].x() == 0
        assert rect1.regionMarkers()[0].y() == 0

        assert rect1.regionMarkers()[0].marker() == 10

    def test_region_markerposition_pos_size(self):
        # scaled version
        rect1 = mt.createRectangle(
            pos=[2.0, -2.75],
            size=[4.0, 2.5],
            isClosed=True,
            marker=5,
            markerPosition=[0.0, 0.0],
        )
        assert rect1.regionMarkers()[0].x() == 0
        assert rect1.regionMarkers()[0].y() == 0
        assert rect1.regionMarkers()[0].marker() == 5


class TestCreateCircle(unittest.TestCase):
    def test_default_create(self):
        circle = mt.createCircle(
            pos=[0.0, 0.0],
            radius=1.0,
            marker=1
        )
        assert abs(circle.regionMarkers()[0].dist(circle.node(0).pos()) - 0.001) < 1e-8
        assert circle.regionMarkers()[0].marker() == 1

    def test_default_create_with_scaling(self):
        circle = mt.createCircle(
            pos=[0.0, 0.0],
            radius=3.0,
            marker=6
        )
        assert abs(circle.regionMarkers()[0].dist(circle.node(0).pos()) - 0.001) < 1e-8
        assert circle.regionMarkers()[0].marker() == 6

    def test_create_with_custom_markerPosition(self):
        circle = mt.createCircle(
            pos=[0.0, 0.0],
            radius=3.0,
            marker=9,
            markerPosition=[2.0, -1.0],
        )
        assert circle.regionMarkers()[0].x() == 2.0
        assert circle.regionMarkers()[0].y() == -1
        assert circle.regionMarkers()[0].marker() == 9


class TestCreatePolygon(unittest.TestCase):
    def test_default_create(self):
        polygon = mt.createPolygon(
            [[-1.0, -1.0], [-1.0, 1.0], [1.0, 1.0], [1.0, -1]],
            isClosed=True,
            marker=1,
        )

        assert abs(polygon.regionMarkers()[0].dist(polygon.node(0).pos()) - 0.001) < 1e-8
        assert polygon.regionMarkers()[0].marker() == 1

    def test_default_create_different_center(self):
        polygon = mt.createPolygon(
            [[0.0, 0.0], [1.0, 0.0], [1.0, -1.0], [0.0, -1]],
            isClosed=True,
            marker=2,
        )
        assert abs(polygon.regionMarkers()[0].dist(polygon.node(0).pos()) - 0.001) < 1e-8
        assert polygon.regionMarkers()[0].marker() == 2

    def test_create_with_custom_markerPosition(self):
        polygon = mt.createPolygon(
            [[0.0, 0.0], [1.0, 0.0], [1.0, -1.0], [0.0, -1]],
            isClosed=True,
            marker=7,
            markerPosition=[0.1, -0.1],
        )

        assert polygon.regionMarkers()[0].x() == 0.1
        assert polygon.regionMarkers()[0].y() == -0.1
        assert polygon.regionMarkers()[0].marker() == 7

class Test3DMerge(unittest.TestCase):
    def test_cubeBasics(self):
        plc = mt.createCube()
        for i, b in enumerate(plc.boundaries()):
            b.setMarker(i+1)

        mesh = mt.createMesh(plc)

        for marker in pg.unique(pg.sort(plc.boundaryMarkers())):
            b1 = plc.boundaries(plc.boundaryMarkers() == marker)[0]
            b2 = mesh.boundaries(mesh.boundaryMarkers() == marker)[0]
            
            np.testing.assert_array_equal(b1.norm(), b2.norm())

    def test_cube_cube_same(self):
        c1 = mt.createCube()
        c2 = mt.createCube()
        m = mt.mergePLC3D([c1, c2])
        self.assertEqual(c1.nodeCount(), m.nodeCount())
        self.assertEqual(c1.boundaryCount(), m.boundaryCount())

    def test_cube_cube_equalface(self):
        w = mt.createCube(marker=1)
        c = mt.createCube(marker=2)
        c.translate([c.xmax()-w.xmin(), 0.0])

        w = mt.mergePLC3D([w, c])
        self.assertEqual(w.nodeCount(), 8+4)
        self.assertEqual(w.boundaryCount(), 6+5)

        c = mt.createCube(marker=3)
        c.translate([0.0, w.ymax()-c.ymin(), 0.0])
        w = mt.mergePLC3D([w, c])
        self.assertEqual(w.nodeCount(), 8+4+4)
        self.assertEqual(w.boundaryCount(), 6+5+5)

        c = mt.createCube(marker=4)
        c.translate([0.0, 0.0, c.zmax()-w.zmin()])
        w = mt.mergePLC3D([c, w])
        self.assertEqual(w.nodeCount(), 8+4+4+4)
        self.assertEqual(w.boundaryCount(), 6+5+5+5)

        c = mt.createCube(marker=5)
        c.translate([0.0, w.ymax()-c.ymin(), c.zmax()-w.zmin()])
        w = mt.mergePLC3D([c, w])
        self.assertEqual(w.nodeCount(), 8+4+4+4+6)
        self.assertEqual(w.boundaryCount(), 6+5+5+5+6)

        c = mt.createCube(marker=6)
        c.translate([0.0, c.ymax()-w.ymin(), c.zmax()-w.zmin()])
        w = mt.mergePLC3D([w, c])
        self.assertEqual(w.nodeCount(), 8+4+4+4+6+0)
        self.assertEqual(w.boundaryCount(), 6+5+5+5+6+3)

        # w.exportPLC('t.poly')
        # pg.show(mt.createMesh(w))

    def test_cube_cube_coplanar_touchface(self):
        w = mt.createCube(marker=1)
        w.scale([2.0, 2.0, 2.0])

        c = mt.createCube(marker=2)
        c.translate([1.5, 0.0, 0.0])
        w = mt.mergePLC3D([w, c])
        self.assertEqual(w.nodeCount(), 8+8)
        self.assertEqual(w.boundaryCount(), 6+5)

        c = mt.createCube(marker=3)
        c.translate([-1.5, 0.0, 0.0])
        w = mt.mergePLC3D([w, c])
        self.assertEqual(w.nodeCount(), 8+8+8)
        self.assertEqual(w.boundaryCount(), 6+5+5)

        c = mt.createCube(marker=4)
        c.translate([0.0, 1.5, 0.0])
        w = mt.mergePLC3D([w, c])
        self.assertEqual(w.nodeCount(), 8+8+8+8)
        self.assertEqual(w.boundaryCount(), 6+5+5+5)

        c = mt.createCube(marker=5)
        c.translate([0.0, 0.0, -1.5])
        w = mt.mergePLC3D([w, c])
        self.assertEqual(w.nodeCount(), 8+8+8+8+8)
        self.assertEqual(w.boundaryCount(), 6+5+5+5+5)

        #pg.show(w)
        # w.exportPLC('t.poly')
        # pg.show(mt.createMesh(w))

    def test_smallcube_in_bigcube(self):
        """
        A small cube in a bigger one, creating two subfaces.
        author: @frodo4fingers
        """
        w = mt.createCube(marker=1)
        c = mt.createCube(size=[0.5, 1.0, 1.0], marker=2)

        w = mt.mergePLC3D([w, c])
        self.assertEqual(w.nodeCount(), 8+8)
        self.assertEqual(w.boundaryCount(), 8)

        # will not work until edge intersection is working
        # d = mt.createCube(size=[0.8, 1.0, 1.0],
        #                   pos=[0.1, 0.0, 1.0],
        #                   marker=3)
        # w = mt.mergePLC3D([w, d])
        # self.assertEqual(w.nodeCount(), 8+8)
        # self.assertEqual(w.boundaryCount(), 8)

        # print(w)
        # pg.show(w)
        # w.exportPLC('w3D_test.w')
        # pg.show(mt.createMesh(w))

    def test_face_in_face(self):
        """Test subface with different marker constructed with hole marker."""
        w = mt.createCube(marker=1, boundaryMarker=1)
        b = w.boundary(2)

        pad = mt.createFacet(mt.createCircle(radius=0.2, segments=12,
                                             isHole=True))
        b2 = pad.boundary(0)

        # rotate to match target norm and pos
        rot = pg.core.getRotation(b2.norm(), b.norm())
        pad.transform(rot)
        pad.translate(b.center())
        # create a boundary with new marker match the hole
        w.copyBoundary(b2)

        w.createBoundary(w.nodes([w.createNode(n.pos()).id() for n in b2.nodes() ]), marker=2)

        # w.exportPLC('pad.poly')
        mesh = mt.createMesh(w)

        np.testing.assert_array_equal(pg.unique(pg.sort(mesh.boundaryMarkers())), [1, 2])

        # print(mesh)
        # mesh.exportBoundaryVTU('b.vtu')
        #pg.show(mesh)

if __name__ == '__main__':
    unittest.main()
