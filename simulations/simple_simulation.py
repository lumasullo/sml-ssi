#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 18 15:58:49 2021

@author: Luciano
"""

import os
from os.path import dirname as up

# use .../sml-ssi main folder as working directory
cwd = os.getcwd()
wdir = up(up(cwd))
os.chdir(wdir)

import numpy as np
import tools.colormaps as cmaps
import tools.tools_simulations as tools
import matplotlib.pyplot as plt

plt.close('all')

#%% Set parameters and initialize arrays

method = 'RASTMIN'
psf_type = 'doughnut'
center_value = False
N = 25000 # detected photons
SBR = 1000 # Signal to Background Ratio
L = 100 # characteristic distance 
K = 25
fwhm = 300 # fwhm of the psf
size_nm = 200 # field of view size (nm)
px_nm = 1 # digital resolution
size = int(size_nm/px_nm)

r0_nm = [1, 1]
samples = 2

extent = [-size_nm/2, size_nm/2, -size_nm/2, size_nm/2]

x = np.arange(-size/2, size/2)
y = np.arange(-size/2, size/2)

Mx, My = np.meshgrid(x, y)
Mr = np.sqrt(Mx**2 + My**2)

pos_nm = tools.ebp_centres(K, L, center=center_value, phi=0, 
                           arr_type='raster scan')

psf = np.zeros((K, size, size)) # array of sequential illuminations

#%% Simulate PSFs

for i in range(K):
    
    psf[i, :, :] = tools.psf(pos_nm[i, :], size_nm, px_nm, fwhm, 
                             psf_type=psf_type)

#%% Simulate n_array, and estimate position through MLE

# "localizations" parameters must match

n_array = tools.sim_exp(psf, r0_nm, N, SBR, size_nm, px_nm, 
                        localizations=samples, DEBUG=False)

r0_mle = tools.mle(n_array, psf, SBR, px_nm=px_nm, prior='rough loc', s=50, 
                   localizations=samples, DEBUG=False)
    
err_array = r0_mle - r0_nm

print('2D error is', np.sqrt((1/2)*np.mean(err_array[:, 0]**2+err_array[:, 1]**2)))

fig, ax = plt.subplots()

# show the first simulated n_array as an example
ax.imshow(n_array[0].reshape(int(np.sqrt(K)), int(np.sqrt(K))), 
          interpolation='None', cmap=cmaps.parula)
    
plt.figure('x estimator')
plt.hist(r0_mle[:, 0])

plt.figure('y estimator')
plt.hist(r0_mle[:, 1])