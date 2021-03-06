#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 13 14:43:20 2021

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
from datetime import datetime
import configparser

#%% Set parameters and initialize arrays

method = 'MINFLUX'
psf_type = 'doughnut'
center_value = True
SBR_array = np.logspace(-0.3, 1.47, num=25) # Signal to Background Ratio
SBR_array = np.append(SBR_array, [5, 10, 100000])
N = 500
L = 150 # distance between beam centers
fwhm = 300 # fwhm of the psf
fov = 0.75*L # fov for the average σ_CRB
size_nm = 1.2*L # field of view size (nm)
step_nm = 1 # digital resolution
size = int(size_nm/step_nm)

K = 4

extent = [-size_nm/2, size_nm/2, -size_nm/2, size_nm/2]

x = np.arange(-size/2, size/2)
y = np.arange(-size/2, size/2)

Mx, My = np.meshgrid(x, y)
Mr = np.sqrt(Mx**2 + My**2)

av_σ_array = np.zeros(len(SBR_array))

pos_nm = tools.ebp_centres(K, L, center=center_value, arr_type='orbit', phi=0)

psf = np.zeros((K, size, size)) # array of sequential illuminations

#%% Simulate PSFs

for i in range(K):
    
    psf[i, :, :] = tools.psf(pos_nm[i, :], size_nm, step_nm, fwhm,
                             psf_type=psf_type)
    
#%% Calculate CRB for different N values

for i, SBR in enumerate(SBR_array):

    σ_CRB, Σ_CRB, Fr, sbr_rel = tools.crb(K, psf, SBR, step_nm, size_nm, N)
    
    mask = tools.create_circular_mask(size, size, center=None, radius=fov/2)
    σ_CRB_cropped = σ_CRB[mask]
    av_σ = np.mean(σ_CRB_cropped)
    av_σ_array[i] = np.mean(σ_CRB_cropped)
    
fig, ax = plt.subplots()

ax.plot(SBR_array[:-3], av_σ_array[:-3])
ax.set_xlabel('SBR')
ax.set_ylabel('average σ_CRB (nm)')


#%% Save results

path = os.getcwd()
filename = r'/minflux_sigma_vs_sbr_L_' + str(L)
folder = r'/results'
np.save(path + folder + filename + '_av_sigma_array.npy', av_σ_array)
np.save(path + folder + filename + '_sbr_array' + '.npy', SBR_array)

SBR = 'variable'

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
'px (nm)': step_nm,
'psf_type': psf_type,
'central excitation': center_value,
'file name': filename}

with open(path + folder + filename + '_params.txt', 'w') as configfile:
    config.write(configfile)