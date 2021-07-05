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

N_array = np.logspace(1.3, 3.477, num=25) # detected photons
SBR = 5 # Signal to Background Ratio
L = 1200 # characteristic distance 
fov = .01*L # fov for the average σ_CRB
K = 64

fwhm = 300 # fwhm of the psf
sigma_psf = fwhm/2.35
px_size_nm = L / (np.sqrt(2)*(np.sqrt(K)-1)) # equivalent px size to RASTMAX

size_nm = 120 # field of view size (nm)
dx_nm = 1 # digital resolution
size = int(size_nm/dx_nm)

extent = [-size_nm/2, size_nm/2, -size_nm/2, size_nm/2]

x = np.arange(-size/2, size/2)
y = np.arange(-size/2, size/2)

Mx, My = np.meshgrid(x, y)
Mr = np.sqrt(Mx**2 + My**2)

#%% Calculate CRB for different N values

av_σ_array = np.zeros(len(N_array))

for i, N in enumerate(N_array):
    
    print(i)

    σ_CRB = tools.crb_camera(K, px_size_nm, sigma_psf, SBR, N, size_nm, dx_nm)

    mask = tools.create_circular_mask(size, size, center=None, radius=fov/2)
    σ_CRB_cropped = σ_CRB[mask]
    av_σ = np.mean(σ_CRB_cropped)
    av_σ_array[i] = av_σ
    
fig, ax = plt.subplots()

ax.plot(N_array, av_σ_array)
ax.set_xlabel('N')
ax.set_ylabel('average σ_CRB (nm)')

#%% Save results

path = os.getcwd()
filename = r'/camera_sigma_vs_n_L_' + str(L)
folder = r'/results'
np.save(path + folder + filename + '_av_sigma_array.npy', av_σ_array)
np.save(path + folder + filename + '_N_array' + '.npy', N_array)

N = 'variable'

config = configparser.ConfigParser()

config['params'] = {

'Date and time': str(datetime.now()),
'N': N,
'SBR': SBR,
'L (nm)': L,
'fov (nm)': fov,
'K': K,
'fwhm (nm)': fwhm,
'size (nm)': size_nm,
'physical px (nm)': px_size_nm,
'simulation px (nm)': dx_nm,
'file name': filename}

with open(path + folder + filename + '_params.txt', 'w') as configfile:
    config.write(configfile)