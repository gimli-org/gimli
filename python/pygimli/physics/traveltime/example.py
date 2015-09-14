#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Some additional infos here?
"""
import os
from pygimli.physics import Refraction


ra = Refraction(os.path.dirname(__file__) + '/example_topo.sgt')
print(ra)
ra.showData()
ra.showVA()
ra.makeMesh()
ra.showMesh()
ra.run(lam=300)
ra.showResult()
