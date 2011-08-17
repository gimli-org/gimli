#! /usr/bin/python
# -*- coding: utf-8 -*-

import sys

sys.path.append( '../' )
import pygimli as g

print g.versionStr()
print "Lin solver autodetection chose: %s " % g.LinSolver( ).solverName()
