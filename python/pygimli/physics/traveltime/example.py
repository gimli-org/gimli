#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Some additional infos here?
"""
import os
from pygimli.physics import Refraction
import matplotlib.pyplot as plt

ra = Refraction(os.path.dirname(__file__) + '/example_topo.sgt')
print(ra)
ra.showVA()

ra.createMesh()
ra.inv.setMaxIter(20)
ra.run(lam=300)

print("model:", ra.model())
print("paraDomain:", ra.paraDomain())

# %%
fig, ax = plt.subplots(nrows=2)
ra.showResult(ax=ax[0], cMin=300, cMax=1500)
ra.showData(ax=ax[1], response=ra.response)

pg.wait()

