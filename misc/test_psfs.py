#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 17 19:06:10 2021

@author: Luciano
"""

import os
wdir = r'/Users/Luciano/Documents/GitHub/ssi-sml'
os.chdir(wdir)

import numpy as np
import tools.colormaps as cmaps
import tools.tools_simulations as tools
import matplotlib.pyplot as plt

plt.close('all')

fwhm = 300 # fwhm of the psf
size_nm = 1200 # field of view size (nm)
step_nm = 1 # digital resolution
size = int(size_nm/step_nm)

pos_nm = np.array([10, 10])

x = np.arange(-size/2, size/2)
y = np.arange(-size/2, size/2)

Mx, My = np.meshgrid(x, y)
Mr = np.sqrt(Mx**2 + My**2)

gaussian_psf = tools.psf(pos_nm, size_nm, step_nm, fwhm, psf_type='gaussian')

extent = [-size_nm/2, size_nm/2, -size_nm/2, size_nm/2]

plt.figure()
plt.imshow(gaussian_psf, extent=extent, interpolation='None', 
           cmap=cmaps.parula)

plt.figure()
plt.plot(gaussian_psf[int(size/2), :])

doughnut_psf = tools.psf(pos_nm, size_nm, step_nm, fwhm, psf_type='doughnut')

extent = [-size_nm/2, size_nm/2, -size_nm/2, size_nm/2]

plt.figure()
plt.imshow(doughnut_psf, extent=extent, interpolation='None', 
           cmap=cmaps.parula)

plt.figure()
plt.plot(doughnut_psf[int(size/2), :])

#%% Load real PSF

from skimage import io

im = io.imread('AVG_PSF.tif')
imarray = np.array(im)
psf_exp = imarray.astype(float)

#%% 2 PSFs

#doughnut_psf0 = tools.psf(np.array([10, 40]), size_nm, step_nm, fwhm, psf_type='doughnut')
#doughnut_psf1 = tools.psf(np.array([80, 10]), size_nm, step_nm, fwhm, psf_type='doughnut')
#
#mixed_psf = doughnut_psf0 + doughnut_psf1 
#
#plt.figure()
#plt.imshow(mixed_psf, extent=extent, interpolation='None', 
#           cmap=cmaps.parula)
#
#plt.figure()
#plt.plot(mixed_psf[int(size/2), :])
