# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 18:40:51 2023

@author: Guenther.T
"""
import numpy as np
import matplotlib.pyplot as plt

a = 1.0
bvec = [0.6, 1.2, 1.8, 2.4, 3.0]
z = np.arange(0, 2., 0.01)
# %%
fig, ax = plt.subplots()
for b in [0.6, 1.2, 1.8, 2.4, 3.0]:
    sab = np.sqrt(a**2+b**2)
    k = np.pi * b / (1 - b/sab)
    # s = (1/np.sqrt(b**2+z**2) - 1/np.sqrt(sab**2+z**2)) * 2 / 2 / np.pi * k
    s = (z/np.sqrt(b**2+z**2)**3 - z/np.sqrt(sab**2+z**2)**3) * 2 / 2 / np.pi * k
    # s = (z/np.sqrt(b**2+z**2)**3/b - z/np.sqrt(sab**2+z**2)**3/sab) * 2 / 2 / np.pi * k
    # s = (1/b-b/np.sqrt(b**2+z**2) - 1/sab + sab/np.sqrt(sab**2+z**2)) * 2 / 2 / np.pi * k
    ax.plot(s, z)
    # ax.plot(s/max(s), z)

ax.invert_yaxis()
ax.grid(True)
