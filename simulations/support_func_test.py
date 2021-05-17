#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 13:16:20 2020

@author: Luciano
"""

import numpy as np
import matplotlib.pyplot as plt

plt.close('all')


L = 100
R = L/2
D = 300

x = np.linspace(-D, D, 1000)
y = np.linspace(-D, D, 1000)

xx, yy = np.meshgrid(x, y)

b = 4
c = np.sqrt(6)

G = np.exp(-(xx**b)/((R*c)**b)-(yy**b)/((R*c)**b))

fig, ax = plt.subplots()


Gfig = ax.imshow(G, interpolation=None, 
                   extent=[-D, D, -D, D])

ax.set_ylabel('y (nm)')
ax.set_xlabel('x (nm)')

cbar = fig.colorbar(Gfig, ax=ax)
cbar.ax.set_ylabel('Support value')

circ = plt.Circle((0,0), radius=L/2, zorder=10, linestyle='--', facecolor='None', edgecolor='k')
ax.add_patch(circ)

plt.figure('Support function 1D cut')
plt.plot(x, G[499, :])
plt.grid()
plt.vlines(-L/2, 0, 1, linestyles='dashed')
plt.vlines(L/2, 0, 1, linestyles='dashed')
plt.xlabel('x [nm]')
plt.ylabel('Support value (AU)')



