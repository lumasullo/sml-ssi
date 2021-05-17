#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 13:54:33 2020

@author: Luciano

This scripts uses the functions in tools.tools_simulations to simulate a
desired number (samples) of p-MINFLUX experiments, either with or without
blinking dynamics.

At the end it saves, plots and summarizes the results of the simulations.

"""

import os 
path = os.getcwd() # get current directory 
wdir = os.path.abspath(os.path.join(path, os.pardir)) # gets parent directory 
os.chdir(wdir)

import numpy as np
import matplotlib.pyplot as plt
import tools.tools_simulations as tools
from datetime import datetime
import configparser
from PIL import Image

DEBUG = True
plt.close('all')

#%% Number of samples of the simulation

samples = 1000

# folder
now = str(datetime.now()).replace('-', '')
now = now.replace(':', ' ')
now = now.replace('.', '_')
root = r'/Users/Luciano/Documents/GitHub/MINFLUX-software/minflux_estimator/Results/'
folder = root + now + ' p_minflux'

os.mkdir(folder)

print(datetime.now(), 'Successfully created the directory {}'.format(folder))

#%%  simulated simplified (not p nor cw) MINFLUX experiment parameters

K = 4   # number of excitation beams

Ns = 90 # photons detected
Nb = 10

SBR = Ns/Nb # Signal to Background Ratio
L = 100 # ditance between beam centers

fov_center = [0, 0] # field of view center
r0_nm = np.array([29, 33])  # position of the emmiter molecule in nm
#r0_nm = L/2 * np.array([np.cos(θ), np.sin(θ)])

#θ = np.pi/4
print('r0_nm', r0_nm)
beamcenter = False

size_nm = 250 # field of view size (nm)
px_nm = 1
size = int(size_nm/px_nm)

# Emitter position in grid coordinates
r0 = tools.spaceToIndex(r0_nm, size_nm, px_nm) 
print('r0', r0)
#r0 = [400, 400]
# EBP centers 
pos_nm = tools.ebp_centres(K, L, center=beamcenter)
#print(pos_nm)

#pos_nm = np.array([[-50, -50], [-50, -17], [-50, 17], [-50, 50],
#                   [-17, -50], [-17, -17], [-17, 17], [-17, 50],
#                   [17, -50], [17, -17], [17, 17], [17, 50],
#                   [50, -50], [50, -17], [50, 17], [50, 50]])

M_p, Tlife, factor, dt, cycle_time, factor = ([None] for i in range(6))
tot_time, exp_time, M = ([None] for i in range(3))

#%% Simulate PSFs

PSFs = np.zeros((K, size, size))
extent = [-size_nm/2, size_nm/2, -size_nm/2, size_nm/2]

for i in range(K):
    PSFs[i, :, :] = tools.psf(pos_nm[i, :], size_nm, px_nm, fov_center, d='donut')
    
    if DEBUG:
        
        filename = folder + '/sim_results' + 'PSF_' + str(i)
        data = PSFs[i, :, :]
        result = Image.fromarray(data)
        result.save(r'{}.tiff'.format(filename))
    
if DEBUG:
    
    pass
#    fig, axes = plt.subplots(2,2)
#    for i, ax in enumerate(axes.flat):
#        ax.set(xlabel='x (nm)', ylabel='y (nm)')
#        ax.imshow(PSFs[i, :, :], interpolation=None, extent=extent)
#    
#    plt.tight_layout()
            
#%% Save simulation parameters
    
filename = folder + '/sim_params'

config = configparser.ConfigParser()

config['simplified_minflux simulation parameters'] = {

    'Date and time': str(datetime.now()),
    'samples': samples,
    'p-minflux equivalent cycles': M_p,
    'Base time resolution (p-minflux cycle time) [ns]': dt,
    'cw-minflux cycle time [ns]': cycle_time,
    'Number of expositions': K,
    'Time per exposition [ns]': exp_time,
    'Fluorescence lifetime of the emitter [ns]': Tlife,
    'Number of cw-minflux cycle': M,
    'N photons (signal)': Ns,
    'N photons (bkg)': Nb,
    'SBR': SBR,
    'L EBP parameter [nm]': L,
    'Beam center': beamcenter,
    'Position of the emitter [nm]': r0_nm,
    'tot_time [μs]': tot_time}

with open(filename + '.txt', 'w') as configfile:
    config.write(configfile)
        
#%% do MINFLUX analysis in loop
    
fail_count = 0
r0_est_nm_array = np.zeros((2, samples))
loglikel = np.zeros((samples, size, size))


for i in range(samples):
    
    print('Sample', i)
    
    # analyze data with p-MINFLUX algorithm
    params = [PSFs, r0, SBR, Ns, Nb, M_p, Tlife, factor, dt, cycle_time]
    n_array = tools.sim_exp('simplified', None, *params)
    print(n_array)
            
    r0_est, loglikel[i] = tools.pos_MINFLUX(n_array, PSFs, SBR=SBR, prior=None, 
                                            L=L, px_nm=px_nm, DEBUG=True)
    r0_est_nm = tools.indexToSpace(r0_est, size_nm, px_nm)
    print('r0_est_nm', r0_est_nm)


    r0_est_nm_array[:, i] = r0_est_nm
        
#%% Evaluate and save results        
        
print(r0_est_nm_array)
r0_est_nm_array = r0_est_nm_array[~np.isnan(r0_est_nm_array)] 
r0_est_nm_array = r0_est_nm_array.reshape(2, samples-fail_count)
   
mean = np.mean(r0_est_nm_array, axis=1)
std = np.std(r0_est_nm_array, axis=1)
print('r0_est_nm mean', mean)
print('r0_est_nm std', std)

x_est_array = r0_est_nm_array[0, :]
y_est_array = r0_est_nm_array[1, :]

cov = np.cov(x_est_array, y_est_array)

x = r0_nm[0]
y = r0_nm[1]

err_x_array = np.abs(x_est_array - x)
err_y_array = np.abs(y_est_array - y)

print('2D error is', np.sqrt((1/2)*np.mean(err_x_array**2+err_y_array**2)))

print('Number of failed simulations', fail_count, 'out of', samples)
print('Failed simulations should be kept below 5%')

plt.figure('Histogram x')

nbins=50
plt.hist(r0_est_nm_array[0, :], bins=nbins)

plt.xlabel('x estimator (nm)')
plt.ylabel('Counts')

plt.figure('Histogram y')

nbins=50
plt.hist(r0_est_nm_array[1, :], bins=nbins)

plt.xlabel('y estimator (nm)')
plt.ylabel('Counts')

plt.figure()

l_mean = np.mean(loglikel, axis=0)
fig, ax = plt.subplots()

lfig = ax.imshow(l_mean, interpolation=None, 
                 extent=[-size_nm/2, size_nm/2, -size_nm/2, size_nm/2])

cbar = fig.colorbar(lfig, ax=ax)
cbar.ax.set_ylabel('log-likelihood')

ax.set_xlabel('x (nm)')
ax.set_ylabel('y (nm)')

circ = plt.Circle((0,0), radius=L/2, zorder=10, 
                  facecolor='None', edgecolor='w')
ax.add_patch(circ)

plt.figure()
xrange = np.linspace(-size_nm/2, size_nm/2, size)
plt.plot(xrange, l_mean[int(size/2), :])


        
filename = folder + '/sim_results' + 'logL_'
data = l_mean
result = Image.fromarray(data)
result.save(r'{}.tiff'.format(filename))

#TODO: save as .tiff

##%% Save simulation results
#    
#filename = folder + '/sim_results'
#
#config = configparser.ConfigParser()
#
#config['Simulation results'] = {
#
#    'Date and time': str(datetime.now()),
#    'err_mean': err_mean,
#    'err_std [nm]': err_std,
#    'err_1d [nm]': err_1d}
#
#with open(filename + '.txt', 'w') as configfile:
#    config.write(configfile)
#    
#np.save(filename + '_r0_est_nm_array', r0_est_nm_array)