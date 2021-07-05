#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 13 14:43:20 2021

@author: Luciano A. Masullo
"""

import os
wdir = r'/Users/Luciano/Documents/GitHub/sml-ssi'
os.chdir(wdir)

import numpy as np
import tools.colormaps as cmaps
import tools.tools_simulations as tools
import matplotlib.pyplot as plt
from datetime import datetime
import configparser

#%% Set parameters and initialize arrays

method = 'CAMERA'

N = 500 # detected photons
SBR = 5 # Signal to Background Ratio
L = 600 # characteristic distance 
fov = .75*L # fov for the average σ_CRB
K = 16

fwhm = 300 # fwhm of the psf
sigma_psf = fwhm/2.35
px_size_nm = L / (np.sqrt(2)*(np.sqrt(K)-1)) # equivalent px size to RASTMAX

size_nm = 700 # field of view size (nm)
dx_nm = 1 # digital resolution
size = int(size_nm/dx_nm)

extent = [-size_nm/2, size_nm/2, -size_nm/2, size_nm/2]

x = np.arange(-size/2, size/2)
y = np.arange(-size/2, size/2)

Mx, My = np.meshgrid(x, y)
Mr = np.sqrt(Mx**2 + My**2)

#%% Calculate CRB and plot

σ_CRB = tools.crb_camera(K, px_size_nm, sigma_psf, SBR, N, size_nm, dx_nm)


fig, ax = plt.subplots()
fig.suptitle(method + ' CRB')

crbfig = ax.imshow(σ_CRB, interpolation=None, 
                   extent=[-size_nm/2, size_nm/2, -size_nm/2, size_nm/2], 
                   cmap=cmaps.parula, vmin=None, vmax=None)

ax.set_ylabel('y (nm)')
ax.set_xlabel('x (nm)')
ax.set_xlim(-size_nm/2, size_nm/2)
ax.set_ylim(-size_nm/2, size_nm/2)

cbar = fig.colorbar(crbfig, ax=ax)
cbar.ax.set_ylabel('$σ_{CRB}$ (nm)')
    
mask = tools.create_circular_mask(size, size, radius=fov/2)
σ_CRB_cropped = σ_CRB[mask]
av_sigma = np.mean(σ_CRB_cropped)

print('Average precision is', np.around(av_sigma, 2), ' nm')

#%% Save results

#path = os.getcwd()
#filename = r'/rastmax_crb'
#folder = r'/Results'
#np.save(path + folder + filename + '_σ_CRB', σ_CRB)
#
#config = configparser.ConfigParser()
#
#config['params'] = {
#
#'Date and time': str(datetime.now()),
#'N': N,
#'SBR': SBR,
#'L (nm)': L,
#'fov (nm)': fov,
#'K': K,
#'fwhm (nm)': fwhm,
#'size (nm)': size_nm,
#'px (nm)': step_nm,
#'psf_type': psf_type,
#'central excitation': center_value,
#'file name': filename}
#
#with open(path + folder + filename + '_params.txt', 'w') as configfile:
#    config.write(configfile)