#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 18 15:58:49 2021

@author: Luciano
"""

import os
wdir = r'/Users/Luciano/Documents/GitHub/ssi-sml'
os.chdir(wdir)

import numpy as np
import tools.colormaps as cmaps
import tools.tools_simulations as tools
import matplotlib.pyplot as plt

#%% Set parameters and initialize arrays

method = 'RASTMIN'
psf_type = 'doughnut'
center_value = False
N = 2000 # detected photons
SBR = 5 # Signal to Background Ratio
L = 100 # characteristic distance 
K = 81
fov = .75*L # fov for the average σ_CRB
fwhm = 300 # fwhm of the psf
size_nm = 200 # field of view size (nm)
px_nm = .2 # digital resolution
size = int(size_nm/px_nm)

r0_nm = [0, 0]
samples = 100

extent = [-size_nm/2, size_nm/2, -size_nm/2, size_nm/2]

x = np.arange(-size/2, size/2)
y = np.arange(-size/2, size/2)

Mx, My = np.meshgrid(x, y)
Mr = np.sqrt(Mx**2 + My**2)

pos_nm = tools.ebp_centres(K, L, center=center_value, phi=0, 
                           arr_type='raster scan')

#pos_nm = tools.ebp_centres(K, L, center=center_value, phi=0, 
#                           arr_type='orbit')

psf = np.zeros((K, size, size)) # array of sequential illuminations

#%% Simulate PSFs

for i in range(K):
    
    psf[i, :, :] = tools.psf(pos_nm[i, :], size_nm, px_nm, fwhm, 
                             psf_type=psf_type)

r0_mle_array = np.zeros((samples, 2))

#for i in range(samples):
#
#    n_array = tools.sim_exp(psf, r0_nm, N, SBR, size_nm, step_nm, 
#                            localizations=1, DEBUG=False)
#    
#    r0_mle = tools.mle(n_array, psf, SBR, px_nm=1, prior='rough loc', s=50, 
#                       localizations=1, DEBUG=False)
#    
#    print(r0_mle)
#    
#    r0_mle_array[i, :] = r0_mle

n_array = tools.sim_exp(psf, r0_nm, N, SBR, size_nm, px_nm, 
                        localizations=100, DEBUG=False)

r0_mle = tools.mle(n_array, psf, SBR, px_nm=px_nm, prior=None, s=None, 
                   localizations=100, DEBUG=False)
    
r0_mle_array = r0_mle

err_array = r0_mle_array - r0_nm

print('2D error is', np.sqrt((1/2)*np.mean(err_array[:, 0]**2+err_array[:, 1]**2)))

fig, ax = plt.subplots()

ax.imshow(n_array[0].reshape(9, 9), interpolation='None', cmap=cmaps.parula)
    
plt.figure('x estimator')
plt.hist(r0_mle_array[:, 0])

plt.figure('y estimator')
plt.hist(r0_mle_array[:, 1])