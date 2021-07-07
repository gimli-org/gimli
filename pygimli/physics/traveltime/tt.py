#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Plotting functions for traveltime."""
# general purpose
import pygimli as pg
from .TravelTimeManager import TravelTimeManager


def simulate(mesh, scheme, slowness=None, vel=None, **kwargs):
    """Simulate traveltime data."""
    mgr = TravelTimeManager()
    return mgr.simulate(mesh, scheme, slowness=slowness, vel=vel, **kwargs)


simulate.__doc__ = TravelTimeManager.__doc__


class DataContainerTT(pg.DataContainer):
    """Data Container for traveltime."""

    def __init__(self, data=None, **kwargs):
        """Initialize empty data container, load or copy existing."""
        if isinstance(data, pg.DataContainer):
            super().__init__(data, **kwargs)
            self.registerSensorIndex("s")
            self.registerSensorIndex("g")
            self.setSensorIndexOnFileFromOne(True)
        else:
            super().__init__(**kwargs)
            self.registerSensorIndex("s")
            self.registerSensorIndex("g")
            if isinstance(data, str):
                self.load(data)
