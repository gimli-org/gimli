"""Tests for pygimli.meshtools.polytools
"""
import pygimli.meshtools as mt


class TestCreateRectangle(object):
    def test_region_marker_position_basics(self):
        rect1 = mt.createRectangle(
            start=[0.0, 0.0],
            end=[2.0, -1.0],
            isClosed=True,
            marker=1,
        )
        # by default the region marker should be located at
        # sPos + (ePos - sPos) * 0.2)
        assert rect1.regionMarker()[0].x() == 0.4
        assert rect1.regionMarker()[0].y() == -0.2
        assert rect1.regionMarker()[0].marker() == 1

    def test_region_marker_position_translation_scale(self):
        rect1 = mt.createRectangle(
            pos=[1.0, 2.0],
            size=[2.0, 4],
            isClosed=True,
            marker=20,
        )
        # by default the region marker should be located at
        # sPos + (ePos - sPos) * 0.2)
        assert rect1.regionMarker()[0].x() == 0.4
        assert rect1.regionMarker()[0].y() == 4 - 4 * 0.2
        assert rect1.regionMarker()[0].marker() == 20

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
        assert rect1.regionMarker()[0].x() == -1.2
        assert rect1.regionMarker()[0].y() == -0.6

        assert rect1.regionMarker()[0].marker() == 1
        assert rect2.regionMarker()[0].marker() == 2

        assert rect1.regionMarker()[0].x() == rect1.regionMarker()[0].x()
        assert rect2.regionMarker()[0].x() == rect2.regionMarker()[0].x()

        assert rect1.regionMarker()[0].y() == rect1.regionMarker()[0].y()
        assert rect2.regionMarker()[0].y() == rect2.regionMarker()[0].y()

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
        assert rect1.regionMarker()[0].x() == 4 * 0.2
        assert rect1.regionMarker()[0].y() == -1.5 - 2.5 * 0.2

        assert rect1.regionMarker()[0].marker() == 1
        assert rect2.regionMarker()[0].marker() == 2

        assert rect1.regionMarker()[0].x() == rect1.regionMarker()[0].x()
        assert rect2.regionMarker()[0].x() == rect2.regionMarker()[0].x()

        assert rect1.regionMarker()[0].y() == rect1.regionMarker()[0].y()
        assert rect2.regionMarker()[0].y() == rect2.regionMarker()[0].y()

    def test_region_markerposition_start_end(self):
        rect1 = mt.createRectangle(
            start=[0.0, -1.5],
            end=[4.0, -4.0],
            isClosed=True,
            marker=10,
            markerPosition=[0.0, 0.0],
        )
        assert rect1.regionMarker()[0].x() == 0
        assert rect1.regionMarker()[0].y() == 0

        assert rect1.regionMarker()[0].marker() == 10

    def test_region_markerposition_pos_size(self):
        # scaled version
        rect1 = mt.createRectangle(
            pos=[2.0, -2.75],
            size=[4.0, 2.5],
            isClosed=True,
            marker=5,
            markerPosition=[0.0, 0.0],
        )
        assert rect1.regionMarker()[0].x() == 0
        assert rect1.regionMarker()[0].y() == 0
        assert rect1.regionMarker()[0].marker() == 5


class TestCreateCircle(object):
    def test_default_create(self):
        circle = mt.createCircle(
            pos=[0.0, 0.0],
            radius=1.0,
            marker=1
        )
        assert circle.regionMarker()[0].x() == 0
        assert circle.regionMarker()[0].y() == 0
        assert circle.regionMarker()[0].marker() == 1

    def test_default_create_with_scaling(self):
        circle = mt.createCircle(
            pos=[0.0, 0.0],
            radius=3.0,
            marker=6
        )
        assert circle.regionMarker()[0].x() == 0
        assert circle.regionMarker()[0].y() == 0
        assert circle.regionMarker()[0].marker() == 6

    def test_create_with_custom_markerPosition(self):
        circle = mt.createCircle(
            pos=[0.0, 0.0],
            radius=3.0,
            marker=9,
            markerPosition=[2.0, -1.0],
        )
        assert circle.regionMarker()[0].x() == 2.0
        assert circle.regionMarker()[0].y() == -1
        assert circle.regionMarker()[0].marker() == 9


class TestCreatePolygon(object):
    def test_default_create(self):
        polygon = mt.createPolygon(
            [[-1.0, -1.0], [-1.0, 1.0], [1.0, 1.0], [1.0, -1]],
            isClosed=True,
            marker=1,
        )

        assert polygon.regionMarker()[0].x() == 0
        assert polygon.regionMarker()[0].y() == 0
        assert polygon.regionMarker()[0].marker() == 1

    def test_default_create_different_center(self):
        polygon = mt.createPolygon(
            [[0.0, 0.0], [1.0, 0.0], [1.0, -1.0], [0.0, -1]],
            isClosed=True,
            marker=2,
        )

        assert polygon.regionMarker()[0].x() == 0.5
        assert polygon.regionMarker()[0].y() == -0.5
        assert polygon.regionMarker()[0].marker() == 2

    def test_create_with_custom_markerPosition(self):
        polygon = mt.createPolygon(
            [[0.0, 0.0], [1.0, 0.0], [1.0, -1.0], [0.0, -1]],
            isClosed=True,
            marker=7,
            markerPosition=[0.1, -0.1],
        )

        assert polygon.regionMarker()[0].x() == 0.1
        assert polygon.regionMarker()[0].y() == -0.1
        assert polygon.regionMarker()[0].marker() == 7
