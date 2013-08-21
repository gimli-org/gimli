#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

Gimli Basics
------------

*This introductory sentence should state the intent and goal of the tutorial. Keep it brief.
This next block should state any assumptions that you the writer are making. Present them in list form.*

This is the first step for the modelling tutorial where we actually want to use gimli. And show some necessary basics.

The modelling as well as the inversion part of gimli is usually based on discretization, so we will handle this in here.  

First, the library have to be imported.
To avoid name clashes with other libraries we suggest to import `pygimli` and alias it to a simple abbreviation 'g', e.g., by using

"""
import pygimli as g

"""
Every part of :ref:`sec:api` is bind to python and can be used with the leading g.

For instance get the version for gimli with:
"""

print(g.versionStr())


import pylab as P

"""
.. image:: PLOT2RST.current_figure
"""
    
P.show()
