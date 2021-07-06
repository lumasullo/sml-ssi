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

method = 'OT_min'
psf_type = 'doughnut'
center_value = False
N = 500 # detected photons
SBR = 5 # Signal to Background Ratio
L = 100 # ditance between beam centers
fov = .75*L # fov for the average σ_CRB
fwhm = 300 # fwhm of the psf
size_nm = 300 # field of view size (nm)
step_nm = 1 # digital resolution
size = int(size_nm/step_nm)

K = 50

extent = [-size_nm/2, size_nm/2, -size_nm/2, size_nm/2]

x = np.arange(-size/2, size/2)
y = np.arange(-size/2, size/2)

Mx, My = np.meshgrid(x, y)
Mr = np.sqrt(Mx**2 + My**2)

pos_nm = tools.ebp_centres(K, L, center=center_value, phi=0, arr_type='orbit')

psf = np.zeros((K, size, size)) # array of sequential illuminations

#%% Simulate PSFs

for i in range(K):
        
    psf[i, :, :] = tools.psf(pos_nm[i, :], size_nm, step_nm, fwhm, 
                             psf_type=psf_type)
    
#%% Calculate CRB and plot

σ_CRB, Σ_CRB, Fr, sbr_rel = tools.crb(K, psf, SBR, step_nm, size_nm, N, prior='rough loc')

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

circ = plt.Circle((0,0), radius=L/2, zorder=10, linestyle='--', facecolor='None', edgecolor='k')
ax.add_patch(circ)

markercolor1 = 'wo'
markersize1 = 10
    
ax.plot(pos_nm[:, 0], pos_nm[:, 1], markercolor1, markersize=markersize1,
        markerfacecolor='k', markeredgewidth=1, markeredgecolor='w')
    
mask = tools.create_circular_mask(size, size, radius=fov/2)
σ_CRB_cropped = σ_CRB[mask]
av_sigma = np.mean(σ_CRB_cropped)

print('Average precision is', np.around(av_sigma, 2), ' nm')

#%% Plot  Σ_CRB

from matplotlib.patches import Ellipse

fig, ax = plt.subplots()
fig.suptitle(method + 'Σ CRB')

ax.set_xlim([-L, L])
ax.set_ylim([-L, L])

circ = plt.Circle((0,0), radius=L/2, zorder=10, linestyle='--', facecolor='None', edgecolor='k')
ax.add_patch(circ)

n = 10
xx = np.arange(-fov/2, fov/2, step=n)
yy = np.arange(-fov/2, fov/2, step=n)

for x in xx:
    for y in yy:
        
        [i, j] = tools.space_to_index([x, y], size_nm, step_nm)
        
        w, h, theta = tools.cov_ellipse(Σ_CRB[:, :, i, j], nsig=1)
        print('w, h', w, h)
        
        ellipse = Ellipse((x, y), w, h, theta, color='k', fill=False)
        ax.add_patch(ellipse)

#%% Save results

path = os.getcwd()
filename = r'/otmin_crb'
folder = r'/Results'
np.save(path + folder + filename + '_σ_CRB', σ_CRB)

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