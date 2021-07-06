#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 13 14:43:20 2021

@author: Luciano A. Masullo
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
N = 500 # detected photons
SBR = 5 # Signal to Background Ratio
L = 50 # distance between beam centers
fwhm = 300 # fwhm of the psf
size_nm = 400 # field of view size (nm)
step_nm = 1 # digital resolution
size = int(size_nm/step_nm)

fov_array = np.linspace(10, 300, num=200)
fov_array = np.append(fov_array, [1, 0.75 * L])

fov = 'variable'

K = 4

extent = [-size_nm/2, size_nm/2, -size_nm/2, size_nm/2]

x = np.arange(-size/2, size/2)
y = np.arange(-size/2, size/2)

Mx, My = np.meshgrid(x, y)
Mr = np.sqrt(Mx**2 + My**2)

σ_CRB_array = np.zeros((len(fov_array), size, size))

pos_nm = tools.ebp_centres(K, L, center=center_value, arr_type='orbit', phi=0)

psf = np.zeros((K, size, size)) # array of sequential illuminations
    
#%% Simulate PSFs

for i in range(K):
    
    psf[i, :, :] = tools.psf(pos_nm[i, :], size_nm, step_nm, fwhm,
                             psf_type=psf_type)
 
#%% Calculate CRB and plot
    
σ_CRB, Σ_CRB, Fr, sbr_rel = tools.crb(K, psf, SBR, step_nm, size_nm, N, prior='rough loc')


fig, ax = plt.subplots()
fig.suptitle('MINFLUX CRB')

crbfig = ax.imshow(σ_CRB, interpolation=None, 
                   extent=[-size_nm/2, size_nm/2, -size_nm/2, size_nm/2], 
                   cmap=cmaps.parula, vmin=None, vmax=None)

        
ax.set_ylabel('y (nm)')
ax.set_xlabel('x (nm)')
ax.set_xlim(-size_nm/2, size_nm/2)
ax.set_ylim(-size_nm/2, size_nm/2)


cbar = fig.colorbar(crbfig, ax=ax)
cbar.ax.set_ylabel('CRB_2D [nm]')

circ = plt.Circle((0,0), radius=L/2, zorder=10, linestyle='--', facecolor='None', edgecolor='k')
ax.add_patch(circ)

markercolor1 = 'wo'
markersize1 = 10
    
ax.plot(pos_nm[:, 0], pos_nm[:, 1], markercolor1, markersize=markersize1,
        markerfacecolor='k', markeredgewidth=1, markeredgecolor='w')

#%% Calculate avarage σ_CRB for different fovs

av_σ_array = np.zeros(len(fov_array))

for i, fov in enumerate(fov_array):
    
    mask = tools.create_circular_mask(size, size, center=None, radius=fov/2)
    σ_CRB_cropped = σ_CRB[mask]
    av_σ = np.mean(σ_CRB_cropped)
    av_σ_array[i] = av_σ
    
fig, ax = plt.subplots()

ax.plot(fov_array[:-3], av_σ_array[:-3])
ax.set_xlabel('fov diameter (nm)')
ax.set_ylabel('average σ_CRB (nm)')

#%% Save results

path = os.getcwd()
filename = r'/minflux_sigma_vs_fov_L_' + str(L)
folder = r'/results'
np.save(path + folder + filename + '_av_sigma_array.npy', av_σ_array)
np.save(path + folder + filename + '_fov_array' + '.npy', fov_array)

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

