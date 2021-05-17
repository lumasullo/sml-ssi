#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 19:18:02 2020

@author: Luciano
"""

#import os 
#path = os.getcwd() # get current directory 
#wdir = os.path.abspath(os.path.join(path, os.pardir)) # gets parent directory

import os
wdir = r'/Users/Luciano/Documents/GitHub/ssi-sml'
os.chdir(wdir)

import numpy as np
import tools.colormaps as cmaps
import tools.tools_simulations as tools
import matplotlib.pyplot as plt
from PIL import Image
from datetime import datetime
import time

plt.close('all')
π = np.pi

savefigs = True
fov_only = False

#%% Parameters


# folder
now = str(datetime.now()).replace('-', '')
now = now.replace(':', ' ')
now = now.replace('.', '_')
root = r'/Users/Luciano/Documents/GitHub/MINFLUX-software/CRB/Results/'
folder = root + now + ' p_minflux'

os.mkdir(folder)

print(datetime.now(), 'Successfully created the directory {}'.format(folder))

Ns = 450 # photons detected
Nb = 50
N = Ns + Nb
SBR = Ns/Nb # Signal to Background Ratio
L = 100 # ditance between beam centers
size_nm = 160 # field of view size (nm)
step_nm = 1
size = int(size_nm/step_nm)

#K_array = np.array([4, 5, 6, 7, 8, 9, 10, 20, 50, 100])
K_array = np.array([80])
#K_array = np.array([16])
σ_CRB_array = np.zeros((len(K_array), size, size))

# EBP centers 

fov_center = [0, 0] # field of view center
extent = [-size_nm/2, size_nm/2, -size_nm/2, size_nm/2]

# FOV mask

x = np.arange(-size/2, size/2)
y = np.arange(-size/2, size/2)

Mx, My = np.meshgrid(x, y)
Mr = np.sqrt(Mx**2 + My**2)


#%% Simulate PSFs

for i, K in enumerate(K_array):

    pos_nm = tools.ebp_centres(K, L, center=False, phi=0)
    
#    pos_nm = np.array([[-50, -50], [-50, 0], [-50, 50],
#                       [0, -50], [0, 0], [0, 50],
#                       [50, -50], [50, 0], [50, 50]])
    
#    pos_nm = np.array([[-75, -75], [-75, -25], [-75, 25], [-75, 75],
#                       [-25, -75], [-25, -25], [-25, 25], [-25, 75],
#                       [25, -75], [25, -25], [25, 25], [25, 75],
#                       [75, -75], [75, -25], [75, 25], [75, 75]])
    
#    pos_nm = np.array([[-50, -50], [-50, -17], [-50, 17], [-50, 50],
#                       [-17, -50], [-17, -17], [-17, 17], [-17, 50],
#                       [17, -50], [17, -17], [17, 17], [17, 50],
#                       [50, -50], [50, -17], [50, 17], [50, 50]])
#    
#    pos_nm = np.array([[-100, -100], [-100, -34], [-100, 34], [-100, 100],
#                       [-34, -100], [-34, -34], [-34, 34], [-34, 100],
#                       [34, -100], [34, -34], [34, 34], [34, 100],
#                       [100, -100], [100, -34], [100, 34], [100, 100]])
        
    PSFs = np.zeros((K, size, size))
        
    if savefigs is True:
        
        figname = folder + '/crb_' + 'K= ' + str(K) + '.pdf'
    
    for j in range(K):
        
        PSFs[j, :, :] = tools.psf(pos_nm[j, :], size_nm, step_nm, fov_center, d='gaussian')
        
        if savefigs:
            
            filename = folder + '/sim_results' + 'PSF_' + str(i)
            data = PSFs[j, :, :]
            result = Image.fromarray(data)
            result.save(r'{}.tiff'.format(filename))
            
#    σ_CRB, Σ_CRB, Fr, a, b, c, logL, I_f, σ_CRB2 = tools.crb_minflux(K, PSFs, SBR, step_nm, size_nm, N, 
#                                                                     method='3')
    
    
#    σ_CRB = tools.crb_minflux(K, PSFs, SBR, step_nm, size_nm, N, method='1')
    σ_CRB, Σ_CRB, Fr = tools.crb_minflux(K, PSFs, SBR, step_nm, size_nm, N, method='2')
    
    if fov_only:
        
        σ_CRB[Mr>L/2] = 0
    
    σ_CRB_array[i, :, :] = σ_CRB
    
    if savefigs:
    
        fig, ax = plt.subplots()
#        fig.suptitle('CRB, K=' + str(K))
#        crbfig = ax.imshow(σ_CRB, interpolation=None, 
#                           extent=[-size_nm/2, size_nm/2, -size_nm/2, size_nm/2], 
#                           cmap=cmaps.parula)
        
        crbfig = ax.imshow(σ_CRB, interpolation=None, 
                           extent=[-size_nm/2, size_nm/2, -size_nm/2, size_nm/2], 
                           cmap=cmaps.parula, vmin=1.5, vmax=10)

        
        ax.set_ylabel('y (nm)')
        ax.set_xlabel('x (nm)')
        ax.set_xlim(-L/2, L/2)
        ax.set_ylim(-L/2, L/2)

        
        cbar = fig.colorbar(crbfig, ax=ax)
        cbar.ax.set_ylabel('CRB_2D [nm]')
        
        circ = plt.Circle((0,0), radius=L/2, zorder=10, linestyle='--', facecolor='None', edgecolor='k')
#        circ = plt.Circle((0,0), radius=100/2, zorder=10, linestyle='--', facecolor='None', edgecolor='k')
        ax.add_patch(circ)
        
    markercolor1 = 'wo'
    markersize1 = 10
        
    ax.plot(pos_nm[:, 0], pos_nm[:, 1], markercolor1, markersize=markersize1,
             markerfacecolor='k', markeredgewidth=1, markeredgecolor='w')
    
    plt.savefig(figname, format='pdf')
    time.sleep(.1)
    plt.close('CRB, K=' + str(K))
 
av_sigma = np.sum(σ_CRB_array, axis=(1,2))/(size**2)

#print(σ_CRB[87, 65])

#plt.figure()
#plt.imshow(logL, interpolation=None, 
#           extent=[-size_nm/2, size_nm/2, -size_nm/2, size_nm/2], 
#           cmap=cmaps.parula)

#plt.figure()
#plt.plot(K_array, av_sigma, '-o')
#plt.xlabel('K')
#plt.ylabel('average σ_2D CRB')



