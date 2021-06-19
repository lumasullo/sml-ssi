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

method = 'RASTMIN'
psf_type = 'doughnut'
center_value = False
N = 2000 # detected photons
SBR = 5 # Signal to Background Ratio
L = 50 # characteristic distance 
K = 36
k = int(np.sqrt(K))
fov = .75*L # fov for the average σ_CRB
fwhm = 300 # fwhm of the psf
size_nm = 200 # field of view size (nm)
step_nm = 5 # digital resolution
size = int(size_nm/step_nm)

r0_nm = [0, 20]
samples = 2

extent = [-size_nm/2, size_nm/2, -size_nm/2, size_nm/2]

pos_nm = tools.ebp_centres(K, L, center=center_value, phi=0, 
                           arr_type='raster scan')

#pos_nm = tools.ebp_centres(K, L, center=center_value, phi=0, 
#                           arr_type='orbit')

psf = np.zeros((K, size, size)) # array of sequential illuminations

#%% Simulate PSFs

for i in range(K):
    
    psf[i, :, :] = tools.psf(pos_nm[i, :], size_nm, step_nm, fwhm, 
                             psf_type=psf_type)
    
λs = np.zeros(K)
λb = np.zeros(K)
λ = np.zeros(K)
λmle = np.zeros(K)

r0 = tools.space_to_index(r0_nm, size_nm, step_nm)

for i in range(K):
    
    λs[i] = N*(SBR/(SBR+1)) * (psf[i, r0[0], r0[1]]/np.sum(psf, axis=0))[r0[0], r0[1]]
    
λb = np.mean(λs)/(SBR/K)    

for i in range(K):
    
    λ[i] = λs[i] + λb/K

r0_mle_array = np.zeros((samples, 2))

for i in range(samples):

    n_array = tools.sim_exp(psf, r0_nm, N, SBR, size_nm, step_nm, DEBUG=False)
    
    r0_mle = tools.mle(n_array, psf, SBR, px_nm=step_nm, prior=None, s=None, DEBUG=False)
    
    print(r0_mle)
    
    r0_mle_array[i, :] = r0_mle
    

err_array = r0_mle_array - r0_nm

print('2D error is', np.sqrt((1/2)*np.mean(err_array[:, 0]**2+err_array[:, 1]**2)))

fig, ax = plt.subplots()

nplot = ax.imshow(n_array.reshape(k, k), interpolation='None', cmap='hot',
                  extent=[-L/2, L/2, -L/2, L/2])

cbar = fig.colorbar(nplot, ax=ax)
cbar.ax.set_ylabel('Counts')

ax.set_xlabel('x (nm)')
ax.set_ylabel('y (nm)')

fig, ax = plt.subplots()

λplot = ax.imshow(λ.reshape(k, k), interpolation='None', cmap='hot',
                  extent=[-L/2, L/2, -L/2, L/2])

cbar = fig.colorbar(λplot, ax=ax)
cbar.ax.set_ylabel('Counts')

ax.set_xlabel('x (nm)')
ax.set_ylabel('y (nm)')

r0 = tools.space_to_index(r0_mle, size_nm, step_nm)

for i in range(K):
    
    λmle[i] = N*(SBR/(SBR+1)) * (psf[i, r0[0], r0[1]]/np.sum(psf, axis=0))[r0[0], r0[1]] + λb/K

fig, ax = plt.subplots()

x = L/np.sqrt(K) * np.arange(np.sqrt(K))

ax.plot(x, n_array.reshape(k, k)[5, :], '-s', label='Data')
ax.plot(x, λ.reshape(k, k)[5, :], label='λ (ground truth)')
ax.plot(x, λmle.reshape(k, k)[5, :], label='$λ_{MLE}$')

ax.legend()

ax.set_xlabel('x (nm)')
ax.set_ylabel('Counts')


#plt.figure('x estimator')
#plt.hist(r0_mle_array[:, 0])
#
#plt.figure('y estimator')
#plt.hist(r0_mle_array[:, 1])