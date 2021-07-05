#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 14 14:38:22 2021

@author: Luciano
"""

import os
wdir = r'/Users/Luciano/Documents/GitHub/sml-ssi'
os.chdir(wdir)

import numpy as np
import tools.colormaps as cmaps
import tools.tools_simulations as tools
import matplotlib.pyplot as plt

size_nm = 200
px_nm = 1

r0_nm = np.array([0, 1])

r0 = tools.space_to_index(r0_nm, size_nm, px_nm)

print(r0)

r0 = np.array([199, 199])

r0_nm = tools.index_to_space(r0, size_nm, px_nm)

print(r0_nm)

r0 = np.array([0, 0])

r0_nm = tools.index_to_space(r0, size_nm, px_nm)

print(r0_nm)

r0 = np.array([0, 199])

r0_nm = tools.index_to_space(r0, size_nm, px_nm)

print(r0_nm)

r0 = np.array([199, 0])

r0_nm = tools.index_to_space(r0, size_nm, px_nm)

print(r0_nm)

psf = tools.psf([0, 0], 3000, 1, 300, 'doughnut')

plt.figure()
plt.imshow(psf)

print(np.sum(psf))
