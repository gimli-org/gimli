#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Plotting functions for traveltime."""
# general purpose
import numpy as np
import pygimli as pg
from .TravelTimeManager import TravelTimeManager

def simulate(mesh, scheme, slowness=None, vel=None, **kwargs):
    """Dummy. To be replaced."""
    mgr = TravelTimeManager()
    return mgr.simulate(mesh, scheme, slowness=slowness, vel=vel, **kwargs)

simulate.__doc__ = TravelTimeManager.__doc__