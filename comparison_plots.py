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
import configparser

plt.close('all')

folder = '/figure_comparison/'

os.chdir(wdir + folder)

σ_CRB = dict()
σ_CRB['minflux'] = np.load('minflux_crb_σ_CRB.npy')
σ_CRB['rastmin'] = np.load('rastmin_crb_σ_CRB.npy')
σ_CRB['ot'] = np.load('ot_crb_σ_CRB.npy')
σ_CRB['minsted'] = np.load('minsted_crb_σ_CRB.npy')
σ_CRB['otmin'] = np.load('otmin_crb_σ_CRB.npy')
σ_CRB['rastmax'] = np.load('rastmax_crb_σ_CRB.npy')

σ_CRB_N = dict()
σ_CRB_N['minflux'] = np.load('minflux_sigma_vs_n_L_50_av_sigma_array.npy')
σ_CRB_N['rastmin'] = np.load('rastmin_sigma_vs_n_L_50_av_sigma_array.npy')
σ_CRB_N['ot'] = np.load('ot_sigma_vs_n_L_300_av_sigma_array.npy')
σ_CRB_N['minsted'] = np.load('ot_sigma_vs_n_L_50_av_sigma_array.npy')
σ_CRB_N['otmin'] = np.load('otmin_sigma_vs_n_L_50_av_sigma_array.npy')
σ_CRB_N['rastmax'] = np.load('rastmax_sigma_vs_n_L_600_av_sigma_array.npy')

N_array = np.load('ot_sigma_vs_n_L_50_N_array.npy')
#N_array1 = np.load('minflux_sigma_vs_n_L_50_N_array.npy')

σ_CRB_sbr = dict()
σ_CRB_sbr['minflux'] = np.load('minflux_sigma_vs_sbr_L_50_av_sigma_array.npy')[:-3]
σ_CRB_sbr['rastmin'] = np.load('rastmin_sigma_vs_sbr_L_50_av_sigma_array.npy')[:-3]
σ_CRB_sbr['ot'] = np.load('ot_sigma_vs_sbr_L_300_av_sigma_array.npy')[:-3]
σ_CRB_sbr['minsted'] = np.load('ot_sigma_vs_sbr_L_50_av_sigma_array.npy')[:-3]
σ_CRB_sbr['otmin'] = np.load('otmin_sigma_vs_sbr_L_50_av_sigma_array.npy')[:-3]
σ_CRB_sbr['rastmax'] = np.load('rastmax_sigma_vs_sbr_L_600_av_sigma_array.npy')[:-3]

sbr_array = np.load('ot_sigma_vs_sbr_L_50_sbr_array.npy')[:-3]
#sbr_array1 = np.load('minflux_sigma_vs_sbr_L_50_sbr_array.npy')
#sbr_array2 = np.load('otmin_sigma_vs_sbr_L_50_sbr_array.npy')

σ_CRB_fov = dict()
σ_CRB_fov['minflux'] = np.load('minflux_sigma_vs_fov_L_50_av_sigma_array.npy')[:-2]
σ_CRB_fov['rastmin'] = np.load('rastmin_sigma_vs_fov_L_50_av_sigma_array.npy')[:-2]
σ_CRB_fov['ot'] = np.load('ot_sigma_vs_fov_L_300_av_sigma_array.npy')[:-2]
σ_CRB_fov['minsted'] = np.load('ot_sigma_vs_fov_L_50_av_sigma_array.npy')[:-2]
σ_CRB_fov['otmin'] = np.load('otmin_sigma_vs_fov_L_50_av_sigma_array.npy')[:-2]
σ_CRB_fov['rastmax'] = np.load('rastmax_sigma_vs_fov_L_600_av_sigma_array.npy')[:-2]

fov_array = np.load('ot_sigma_vs_fov_L_50_fov_array.npy')[:-2]

fig, ax = plt.subplots(2, 2)

#%% Plot 1D cut σ_crb vs x

size_nm = 300
px_nm = 1
y0 = int((size_nm/px_nm)/2)
x = np.arange(-size_nm/2, size_nm/2)

size_nm = 700
x_rastmax = np.arange(-size_nm/2, size_nm/2)

ax[0, 0].plot(x, σ_CRB['minflux'][y0, :], label='MINFLUX') 
ax[0, 0].plot(x, σ_CRB['rastmin'][y0, :], label='RASTMIN') 
ax[0, 0].plot(x, σ_CRB['ot'][y0, :], label='OT') 
ax[0, 0].plot(x, σ_CRB['otmin'][y0, :], label='OTMIN') 
ax[0, 0].plot(x, σ_CRB['minsted'][y0, :], label='MINSTED') 
ax[0, 0].plot(x_rastmax, σ_CRB['rastmax'][y0, :], label='RASTMAX') 

ax[0, 0].set_xlabel('x (nm)')
ax[0, 0].set_ylabel('$σ_{CRB}$ (nm)')

ax[0, 0].set_xscale('linear')
ax[0, 0].set_yscale('linear')

ax[0, 0].set_xlim([-50, 50])
ax[0, 0].set_ylim([0, 10])

#ax_inset = ax[0, 0].inset_axes([0.6, 0.6, 0.5, 0.5])
#ax_inset.plot(x, σ_CRB['minflux'][y0, :], label='MINFLUX') 
#ax_inset.plot(x, σ_CRB['rastmin'][y0, :], label='RASTMIN') 
#ax_inset.plot(x, σ_CRB['ot'][y0, :], label='OT') 
#ax_inset.plot(x, σ_CRB['otmin'][y0, :], label='OTMIN') 
#ax_inset.plot(x, σ_CRB['minsted'][y0, :], label='MINSTED') 
#ax_inset.plot(x_csmlm, σ_CRB['csmlm'][y0, :], label='Confocal-SMLM') 
#
#ax_inset.tick_params(left=False, bottom=False, labelleft=False,
#                     labelbottom=False)
#
#ax_inset.set_xlim([-50, 50])
#ax_inset.set_ylim([1, 10])

plt.tight_layout()

#%% Plot 1D σ vs fov

ax[0, 1].plot(fov_array, σ_CRB_fov['minflux'], label='MINFLUX') 
ax[0, 1].plot(fov_array, σ_CRB_fov['rastmin'], label='RASTMIN') 
ax[0, 1].plot(fov_array, σ_CRB_fov['ot'], label='OT') 
ax[0, 1].plot(fov_array, σ_CRB_fov['otmin'], label='OTMIN') 
ax[0, 1].plot(fov_array, σ_CRB_fov['minsted'], label='MINSTED') 
ax[0, 1].plot(fov_array, σ_CRB_fov['rastmax'], label='RASTMAX') 

ax[0, 1].set_xlabel('FOV')
ax[0, 1].set_ylabel('$<σ_{CRB}>$ (nm)')

ax[0, 1].set_xscale('log')
ax[0, 1].set_yscale('log')

ax[0, 1].set_ylim([0.6, 120])

ax[0, 1].legend(fontsize=8, ncol=2, loc='upper left')

plt.tight_layout()


#%% Plot 1D σ vs sbr

ax[1, 0].plot(sbr_array, σ_CRB_sbr['minflux'], label='MINFLUX') 
ax[1, 0].plot(sbr_array, σ_CRB_sbr['rastmin'], label='RASTMIN') 
ax[1, 0].plot(sbr_array, σ_CRB_sbr['ot'], label='OT') 
ax[1, 0].plot(sbr_array, σ_CRB_sbr['otmin'], label='OTMIN') 
ax[1, 0].plot(sbr_array, σ_CRB_sbr['minsted'], label='MINSTED') 
ax[1, 0].plot(sbr_array, σ_CRB_sbr['rastmax'], label='Confocal-SMLM') 

ax[1, 0].set_xlabel('SBR')
ax[1, 0].set_ylabel('$<σ_{CRB}>$ (nm)')

ax[1, 0].set_xscale('log')
ax[1, 0].set_yscale('log')

#ax[1, 1].legend(fontsize=7)

plt.tight_layout()

#%% Plot 1D σ vs N

ax[1, 1].plot(N_array, σ_CRB_N['minflux'], label='MINFLUX') 
ax[1, 1].plot(N_array, σ_CRB_N['rastmin'], label='RASTMIN') 
ax[1, 1].plot(N_array, σ_CRB_N['ot'], label='OT') 
ax[1, 1].plot(N_array, σ_CRB_N['otmin'], label='OTMIN') 
ax[1, 1].plot(N_array, σ_CRB_N['minsted'], label='MINSTED') 
ax[1, 1].plot(N_array, σ_CRB_N['rastmax'], label='Confocal-SMLM') 

ax[1, 1].set_xlabel('N')
ax[1, 1].set_ylabel('$<σ_{CRB}>$ (nm)')

ax[1, 1].set_xscale('log')
ax[1, 1].set_yscale('log')

#ax[0, 1].legend(fontsize=8)
#
#import string
#
#for n, a in enumerate(ax.flat):
#
#    a.text(-0.2, 1.1, string.ascii_uppercase[n], transform=a.transAxes, 
#            size=10, weight='bold')

plt.tight_layout()
