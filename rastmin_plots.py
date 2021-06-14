#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 15 13:18:38 2021

@author: Luciano
"""

import os
wdir = r'/Users/Luciano/Documents/GitHub/ssi-sml'
os.chdir(wdir)

import numpy as np
import tools.colormaps as cmaps
import tools.tools_simulations as tools
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import configparser

folder = '/figure_rastmin/'

os.chdir(wdir + folder)

N_array = np.load('rastmin_sigma_vs_n_L_50_N_array.npy')

σ_CRB_N = dict()
σ_CRB_N['L50'] = np.load('rastmin_sigma_vs_n_L_50_av_sigma_array.npy')
σ_CRB_N['L100'] = np.load('rastmin_sigma_vs_n_L_100_av_sigma_array.npy')
#σ_CRB_N['L120'] =
σ_CRB_N['L150'] = np.load('rastmin_sigma_vs_n_L_150_av_sigma_array.npy')

sbr_array = np.load('rastmin_sigma_vs_sbr_L_50_sbr_array.npy')

σ_CRB_sbr = dict()
σ_CRB_sbr['L50'] = np.load('rastmin_sigma_vs_sbr_L_50_av_sigma_array.npy')
σ_CRB_sbr['L100'] = np.load('rastmin_sigma_vs_sbr_L_100_av_sigma_array.npy')
#σ_CRB_sbr['L120'] =
σ_CRB_sbr['L150'] = np.load('rastmin_sigma_vs_sbr_L_150_av_sigma_array.npy')

fov_array = np.load('rastmin_sigma_vs_fov_L_50_fov_array.npy')

σ_CRB_fov = dict()
σ_CRB_fov['L50'] = np.load('rastmin_sigma_vs_fov_L_50_av_sigma_array.npy')
σ_CRB_fov['L100'] = np.load('rastmin_sigma_vs_fov_L_100_av_sigma_array.npy')
#σ_CRB_sbr['L120'] =
σ_CRB_fov['L150'] = np.load('rastmin_sigma_vs_fov_L_150_av_sigma_array.npy')

params_file = 'rastmin_crb_params.txt'

crb_map = np.load('rastmin_crb_σ_CRB.npy')

fig, ax = plt.subplots(2, 2)

fig.set_size_inches(6, 4.5)

#%% Plot CRB 2D map

config = configparser.ConfigParser()
config.read(params_file)

L = int(config['params']['L (nm)'])
K = int(config['params']['K'])
center_value = config['params']['central excitation']

if center_value == 'False':
    
    center_value = False

elif center_value == 'True':
    
    center_value = True
    
else:
    
    print('center_value must be True/False')
    
pos_nm = tools.ebp_centres(K, L, center=center_value, phi=0, arr_type='raster scan')

size_nm = 300
vmin = 1.5
vmax = 15

crbplot = ax[0,0].imshow(crb_map, interpolation=None, extent=[-size_nm/2, size_nm/2, -size_nm/2, size_nm/2], 
                         cmap=cmaps.parula, vmin=vmin, vmax=vmax)


ax[0,0].set_ylabel('y (nm)')
ax[0,0].set_xlabel('x (nm)')
ax[0,0].set_xlim(-size_nm/4, size_nm/4)
ax[0,0].set_ylim(-size_nm/4, size_nm/4)

cbar = fig.colorbar(crbplot, ax=ax[0,0])
cbar.ax.set_ylabel('$σ_{CRB}$ (nm)')

cbar.set_ticks([5, 10, 15])

circ = plt.Circle((0,0), radius=L/2, zorder=10, linestyle='--', facecolor='None', edgecolor='k')
ax[0,0].add_patch(circ)

markercolor1 = 'wo'
markersize1 = 5
    
ax[0,0].plot(pos_nm[:, 0], pos_nm[:, 1], markercolor1, markersize=markersize1,
        markerfacecolor='k', markeredgewidth=.5, markeredgecolor='w')

#%% Plot 1D σ vs fov

ax[0, 1].plot(fov_array, σ_CRB_fov['L50'], label='L = 50 nm') 
ax[0, 1].plot(fov_array, σ_CRB_fov['L100'], label='L = 100 nm') 
ax[0, 1].plot(fov_array, σ_CRB_fov['L150'], label='L = 150 nm') 

ax[0, 1].set_xlabel('FOV (nm)')
ax[0, 1].set_ylabel('$<σ_{CRB}>$ (nm)')

ax[0, 1].set_xscale('linear')
ax[0, 1].set_yscale('linear')

ax[0, 1].legend(fontsize=9)

#%% Plot 1D σ vs SBR

ax[1, 0].plot(sbr_array, σ_CRB_sbr['L50']) 
ax[1, 0].plot(sbr_array, σ_CRB_sbr['L100']) 
ax[1, 0].plot(sbr_array, σ_CRB_sbr['L150']) 

ax[1, 0].set_xlabel('SBR')
ax[1, 0].set_ylabel('$<σ_{CRB}>$ (nm)')

ax[1, 0].set_xscale('log')
ax[1, 0].set_yscale('log')

ax[1, 0].set_ylim([0.9, 11])


#for axis in [ax[1, 0].xaxis, ax[1, 0].yaxis]:
#for axis in [ax[1, 0].yaxis]:
#    axis.set_major_formatter(ScalarFormatter())


#ax[1, 0].set_xticks([], minor=True)
#ax[1, 0].set_xticks([1, 2, 3, 5, 10, 20])
#ax[1, 0].set_yticks([], minor=True)
#ax[1, 0].set_yticks([1, 2, 3, 5, 7, 10])


#%% Plot 1D σ vs N

ax[1, 1].plot(N_array, σ_CRB_N['L50']) 
ax[1, 1].plot(N_array, σ_CRB_N['L100']) 
ax[1, 1].plot(N_array, σ_CRB_N['L150']) 

ax[1, 1].set_xlabel('N')
ax[1, 1].set_ylabel('$<σ_{CRB}>$ (nm)')

ax[1, 1].set_xscale('log')
ax[1, 1].set_yscale('log')

#for axis in [ax[1, 1].xaxis, ax[1, 1].yaxis]:
#for axis in [ax[1, 1].yaxis]:
#    axis.set_major_formatter(ScalarFormatter())
    
#ax[1, 1].set_yticks([], minor=True)
#ax[1, 1].set_yticks([1, 2, 3, 5, 7, 10])
#ax[1, 1].set_yticks([1, 2, 4, 7, 14])


#ax[1, 0].set_xticks([], minor=True)
#ax[1, 1].set_xticks([20, 100, 1000])



#%% add subfigures labels

#import string
#
#for n, a in enumerate(ax.flat):
#
#    a.text(-0.1, 1.1, string.ascii_uppercase[n], transform=a.transAxes, 
#            size=10, weight='bold')

plt.tight_layout()

