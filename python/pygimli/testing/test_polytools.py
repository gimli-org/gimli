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

    def test_region_marker_position_translation_scale(self):
        rect1 = mt.createRectangle(
            pos=[1.0, 2.0],
            size=[2.0, 4],
            isClosed=True,
            marker=1,
        )
        # by default the region marker should be located at
        # sPos + (ePos - sPos) * 0.2)
        assert rect1.regionMarker()[0].x() == 0.4
        assert rect1.regionMarker()[0].y() == 4 - 4 * 0.2

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

        assert rect1.regionMarker()[0].x() == rect1.regionMarker()[0].x()
        assert rect2.regionMarker()[0].x() == rect2.regionMarker()[0].x()

        assert rect1.regionMarker()[0].y() == rect1.regionMarker()[0].y()
        assert rect2.regionMarker()[0].y() == rect2.regionMarker()[0].y()

    def test_region_markerposition_start_end(self):
        rect1 = mt.createRectangle(
            start=[0.0, -1.5],
            end=[4.0, -4.0],
            isClosed=True,
            marker=1,
            markerPosition=[0.0, 0.0],
        )
        assert rect1.regionMarker()[0].x() == 0
        assert rect1.regionMarker()[0].y() == 0

    def test_region_markerposition_pos_size(self):
        # scaled version
        rect1 = mt.createRectangle(
            pos=[2.0, -2.75],
            size=[4.0, 2.5],
            isClosed=True,
            marker=2,
            markerPosition=[0.0, 0.0],
        )
        assert rect1.regionMarker()[0].x() == 0
        assert rect1.regionMarker()[0].y() == 0
