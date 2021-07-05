#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 27 16:33:23 2021

@author: Luciano
"""

import os
wdir = r'/Users/Luciano/Documents/GitHub/sml-ssi'
os.chdir(wdir)

import numpy as np
import matplotlib.pyplot as plt

folder = '/figure_camera/'

os.chdir(wdir + folder)

N_array = np.load('rastmax_sigma_vs_n_L_1200_N_array.npy')

σ_CRB_N = dict()
σ_CRB_N['rastmax'] = np.load('rastmax_sigma_vs_n_L_1200_av_sigma_array.npy')
σ_CRB_N['camera'] = np.load('camera_sigma_vs_n_L_1200_av_sigma_array.npy')

#%% Plot 1D σ vs N

fig, ax = plt.subplots(figsize=(4, 4))

ax.plot(N_array, σ_CRB_N['rastmax'], color=u'#4059AD', label='RASTMAX')
ax.plot(N_array, σ_CRB_N['camera'], color=u'#F4B942', label='camera')

ax.set_xlabel('N')
ax.set_ylabel('$σ_{CRB}$ (nm)')

ax.set_xscale('log')
ax.set_yscale('log')

ax.set_ylim([0.9, 50])
ax.set_xlim([N_array[0], N_array[-1]])

plt.legend()

plt.tight_layout()