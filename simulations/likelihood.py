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

plt.close('all')

#%% Set parameters and initialize arrays

method = 'OTMIN'
psf_type = 'doughnut'
center_value = False
N = 500 # detected photons
SBR = 5 # Signal to Background Ratio
L = 100 # characteristic distance 
K = 50
fwhm = 300 # fwhm of the psf
size_nm = 500 # field of view size (nm)
px_nm = 1 # digital resolution
size = int(size_nm/px_nm)

r_nm_0 = [0, 0]
r_nm_1 = [53, 0]

samples = 1000

extent = [-size_nm/2, size_nm/2, -size_nm/2, size_nm/2]

x = np.arange(-size/2, size/2)
y = np.arange(-size/2, size/2)

Mx, My = np.meshgrid(x, y)
Mr = np.sqrt(Mx**2 + My**2)

pos_nm = tools.ebp_centres(K, L, center=center_value, phi=0, arr_type='orbit')
#pos_nm = tools.ebp_centres(K, L, center=center_value, phi=0, arr_type='raster scan')


psf = np.zeros((K, size, size)) # array of sequential illuminations

#%% Simulate PSFs

for i in range(K):
    
    psf[i, :, :] = tools.psf(pos_nm[i, :], size_nm, px_nm, fwhm, 
                             psf_type=psf_type)
    
#%% Simulation for position 1
    
# "localizations" parameters must match

n_array = tools.sim_exp(psf, r_nm_0, N, SBR, size_nm, px_nm, 
                        localizations=samples, DEBUG=False)

r_mle, _, likelihood = tools.mle(n_array, psf, SBR, px_nm=px_nm, 
                                  prior='rough loc', s=50, 
                                  localizations=samples, DEBUG=True)
    
err_array = r_mle - np.array(r_nm_0)

print('2D error is', np.sqrt((1/2)*np.mean(err_array[:, 0]**2+err_array[:, 1]**2)))

av_likelihood = np.mean(likelihood, axis=0)/N

fig, ax = plt.subplots(2, 3)

w = 1
n = int(np.ceil((r_mle[:, 0].max() - r_mle[:, 0].min())/w))

ax[0, 0].hist(r_mle[:, 0], bins=n)
ax[0, 0].set_xlabel('x (nm)')
ax[0, 0].set_ylabel('Counts')

w = 1
n = int(np.ceil((r_mle[:, 1].max() - r_mle[:, 1].min())/w))

ax[0, 1].hist(r_mle[:, 1], bins=n)
ax[0, 1].set_xlabel('y (nm)')
ax[0, 1].set_ylabel('Counts')

lfig = ax[0, 2].imshow(av_likelihood, interpolation='None', extent=extent, 
                       vmin=-4.112, vmax=-3.912)

ax[0, 2].set_xlim(-60, 60)
ax[0, 2].set_ylim(-60, 60)

cbar = fig.colorbar(lfig, ax=ax[0, 2])
cbar.ax.set_ylabel('A.U.')

circ = plt.Circle((0,0), radius=L/2, zorder=10, linestyle='--', facecolor='None', edgecolor='w')
ax[0, 2].add_patch(circ)

xmean_r0 = np.mean(r_mle[:, 0])
ymean_r0 = np.mean(r_mle[:, 1])

#%% Simulation for position 2

n_array = tools.sim_exp(psf, r_nm_1, N, SBR, size_nm, px_nm, 
                        localizations=samples, DEBUG=False)

r_mle, _, likelihood = tools.mle(n_array, psf, SBR, px_nm=px_nm, 
                                  prior='rough loc', s=50, 
                                  localizations=samples, DEBUG=True)
    
err_array = r_mle - np.array(r_nm_0)

print('2D error is', np.sqrt((1/2)*np.mean(err_array[:, 0]**2+err_array[:, 1]**2)))

av_likelihood = np.mean(likelihood, axis=0)/N

w = 3
n = int(np.ceil((r_mle[:, 0].max() - r_mle[:, 0].min())/w))

ax[1, 0].hist(r_mle[:, 0], bins=n)
ax[1, 0].set_xlabel('x (nm)')
ax[1, 0].set_ylabel('Counts')

w = 1
n = int(np.ceil((r_mle[:, 1].max() - r_mle[:, 1].min())/w))

ax[1, 1].hist(r_mle[:, 1], bins=n)
ax[1, 1].set_xlabel('y (nm)')
ax[1, 1].set_ylabel('Counts')

lfig = ax[1, 2].imshow(av_likelihood, interpolation='None', extent=extent,
                       vmin=-3.934, vmax=-3.734)

ax[1, 2].set_xlim(-60, 60)
ax[1, 2].set_ylim(-60, 60)

cbar = fig.colorbar(lfig, ax=ax[1, 2])
cbar.ax.set_ylabel('A.U.')

circ = plt.Circle((0,0), radius=L/2, zorder=10, linestyle='--', facecolor='None', edgecolor='w')
ax[1, 2].add_patch(circ)

xmean_r1 = np.mean(r_mle[:, 0])
ymean_r1 = np.mean(r_mle[:, 1])

plt.tight_layout()

print(xmean_r0)
print(ymean_r0)
print(xmean_r1)
print(ymean_r0)


    