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

#%% Set parameters and initialize arrays

method = 'CAMERA'
#psf_type = 'gaussian'
#center_value = True
N = 500 # detected photons
SBR = 5 # Signal to Background Ratio
L = 1200 # characteristic distance 
K = 64
#fov = .75*L # fov for the average σ_CRB
fwhm = 300 # fwhm of the psf
size_nm = 100 # field of view size (nm)
step_nm = 1 # digital resolution
size = int(size_nm/step_nm)

px_size_nm = L / (np.sqrt(2)*(np.sqrt(K)-1))
sigma_psf = fwhm/2.35

r0_nm = [0, 0]
samples = 1

r0_mle_array = np.zeros((samples, 2))

for i in range(samples):

    n_array = tools.sim_exp_camera(r0_nm, K, px_size_nm, sigma_psf, SBR, N, size_nm, step_nm)
    
    r0_mle = tools.mle_camera(n_array, K, px_size_nm, sigma_psf, SBR, N, size_nm, step_nm)
    
    print(r0_mle)
    
    r0_mle_array[i, :] = r0_mle
    
    
print('mean r0_mle', np.mean(r0_mle_array, axis=0))

err_array = r0_mle_array - r0_nm

print('2D error is', np.sqrt((1/2)*np.mean(err_array[:, 0]**2+err_array[:, 1]**2)))

fig, ax = plt.subplots(figsize=(3,2))

total_size = px_size_nm * np.sqrt(K)
nplot = ax.imshow(n_array.reshape(8, 8), interpolation='None', cmap='gray',
                  extent=[-total_size/2, total_size/2, -total_size/2, total_size/2],
                  vmin=0, vmax=60)

cbar = fig.colorbar(nplot, ax=ax)
cbar.ax.set_ylabel('Counts')

ax.set_xlabel('x (nm)')
ax.set_ylabel('y (nm)')

plt.tight_layout()


    
    