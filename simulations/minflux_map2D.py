#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 15:21:36 2020

@author: Luciano A. Masullo


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
import time

DEBUG = False
saverelTimes = False
plt.close('all')

π = np.pi

# folder
now = str(datetime.now()).replace('-', '')
now = now.replace(':', ' ')
now = now.replace('.', '_')
root = r'/Users/Luciano/Documents/GitHub/MINFLUX-software/orbital_tracking_sim/Results/'
folder = root + now + ''

os.mkdir(folder)

print(datetime.now(), 'Successfully created the directory {}'.format(folder))

savefigs = False

samples = 10

#%% Create 2D grid in (x, y) and (r, φ)

#  fixed parameters
#w0 = 300 # (nm) gaussian beam waist

[xmin, xmax] = [-75, 75]
[ymin, ymax] = [-75, 75]
n = 40 # number of pixels in the map
x_array = np.linspace(xmin, xmax, n)
y_array = np.linspace(ymin, ymax, n)

Mx, My = np.meshgrid(x_array, -y_array)

#%%  simulated p-MINFLUX experiment parameters 

## This parameters are not used if the simulation method is 'simplified'

ti = time.time()

#M_p = int(5*10**7)
M_p = int(2*10**5)  # p-minflux equivalent illumination cicles
dt = 25 # p-minflux cycle time [ns]
cycle_time = 'None' # [ns] cw-minflux cycle time, i.e. time it takes for 4 exposures
K = 4   # number of excitation beams
exptime = 'None'
Tlife = 0.001 # in ns

M = 'None' # number of cycles
ct = 'None' # time points in units of dt = 25 ns

Ns = 90 # photons detected
Nb = 10

SBR = Ns/Nb # Signal to Background Ratio
L = 100 # ditance between beam centers
A = L/2 # equivalent OT radius

fov_center = [0, 0] # field of view center
#r0_nm = np.array([3.35, 3.35])  # position of the emmiter molecule in nm

tot_time = (M_p * dt)/10**3 # [μs]
print('Experiment length', tot_time, 'μs')

size_nm = 1250 # field of view size (nm)
step_nm = 1
size = int(size_nm/step_nm)

τ = np.arange(0, K)/K * dt 
a=0.0 # [ns]
b=dt/K # [ns]    

central = True

# Emitter position in grid coordinates
#r0 = tools.spaceToIndex(r0_nm, size_nm, step_nm) 

# EBP centers 
#pos_nm = tools.beams(K, L, center=central, d='donut') # TODO: clean-up beams function
pos_nm = tools.ebp_centres(K, L, center=central, phi=0)

#%% Simulate PSFs

PSFs = np.zeros((K, size, size))
extent = [-size_nm/2, size_nm/2, -size_nm/2, size_nm/2]

for i in range(K):
    PSFs[i, :, :] = tools.psf(pos_nm[i, :], size_nm, step_nm, fov_center, d='donut')
    
    if saverelTimes:
        
        filename = folder + '/sim_results'
        np.save(filename + 'PSF_' + str(i), PSFs[i, :, :])
    
if DEBUG:
    
    fig, axes = plt.subplots(2,2)
    for i, ax in enumerate(axes.flat):
        ax.set(xlabel='x (nm)', ylabel='y (nm)')
        ax.imshow(PSFs[i, :, :], interpolation=None, extent=extent)
    
    plt.tight_layout()

#%% Save simulation parameters
    
filename = folder + '/sim_params'

config = configparser.ConfigParser()

config['cw_minflux simulation parameters'] = {

    'Date and time': str(datetime.now()),
    'samples': samples,
    'p-minflux equivalent cycles': M_p,
    'Base time resolution (p-minflux cycle time) [ns]': dt,
    'cw-minflux cycle time [ns]': cycle_time,
    'Number of expositions': K,
    'Time per exposition [ns]': exptime,
    'Fluorescence lifetime of the emitter [ns]': Tlife,
    'Number of cw-minflux cycle': M,
    'N photons (signal)': Ns,
    'N photons (bkg)': Nb,
    'SBR': SBR,
    'L EBP parameter [nm]': L,
    'tot_time [μs]': tot_time}

with open(filename + '.txt', 'w') as configfile:
    config.write(configfile)
    
#%% Simulation loop
    
fail_count = 0
err_map_xy = np.zeros((len(x_array), len(y_array)))

for i, x in enumerate(x_array):
    
    for j, y in enumerate(y_array): 
        
        x_est_array = np.zeros(samples)
        y_est_array = np.zeros(samples)
        
        r0_nm = np.array([x, y])  # position of the emmiter molecule in nm
        r0 = tools.spaceToIndex(r0_nm, size_nm, step_nm) 
    
        factor = 1.05
        params = [PSFs, r0, SBR, Ns, Nb, M_p, Tlife, factor, dt, cycle_time]
                    
        for k in range(samples):
            
            print([x, y], 'run', k, 'out of', samples)
            n_array = tools.sim_exp('simplified', None, *params)
            
            r0_est = tools.pos_MINFLUX(n_array, PSFs, SBR=SBR, prior=None)
            r0_est_nm = tools.indexToSpace(r0_est, size_nm, step_nm)
                        
            Ntot = np.sum(n_array)
                    
            x_est_array[k] = r0_est_nm[0]
            y_est_array[k] = r0_est_nm[1]

        err_x_array = np.abs(x_est_array - x)
        err_y_array = np.abs(y_est_array - y)
        
        err_2D = np.sqrt((.5)*np.mean(err_x_array**2 + err_y_array**2, axis=0))
        
        print('err_2d', err_2D)

        err_map_xy[j, i] = err_2D # fix to fill the matriz with (x, y) as cartesian coordinates
    
#%% Plot expositions and emitter positions

fig1, ax1 = plt.subplots()
fig1.set_size_inches(7, 7)

#ax1.scatter(pos_nm[:,0], pos_nm[:,1], color='w', label='Expositions')
ax1.scatter(0, 0, facecolors='none', edgecolors='k')
ax1.scatter(pos_nm[:, 0], pos_nm[:, 1], facecolors='none', edgecolors='w')

ax1.grid(False)

circ = plt.Circle((0,0), radius=A, zorder=10, facecolor='None', edgecolor='w')
ax1.add_patch(circ)

ax1.set_xlabel('x (nm)')
ax1.set_ylabel('y (nm)')

ax1.set_xlim([-A*1.1, A*1.1])
ax1.set_ylim([-A*1.1, A*1.1])

err_map_fig = ax1.imshow(err_map_xy, extent=[xmin, xmax, ymin, ymax])
                               
cbar = fig1.colorbar(err_map_fig, ax=ax1)
cbar.ax.set_ylabel('err_2D [nm]')

plt.tight_layout()

tf = time.time()

simtime = tf-ti

#%% Save simulation results
    
filename = folder + '/sim_results'

config = configparser.ConfigParser()

config['Simulation parameters'] = {

    'Date and time': str(datetime.now()),
    'n (map pixels)': n,
    'xmin (nm)': xmin,
    'xmax (nm)': xmax,
    'ymin (nm)': ymin,
    'ymax (nm)': ymax,
    'N (photons)': Ns + Nb,
    'K': K,
    'SBR': SBR,
    'L [nm]': L,
    'fail count': fail_count,
    'Sim time [s]': simtime,
    'Central doughnut': str(central)
    }

with open(filename + '.txt', 'w') as configfile:
    config.write(configfile)
    
np.save(filename + 'err_map_xy', err_map_xy)





