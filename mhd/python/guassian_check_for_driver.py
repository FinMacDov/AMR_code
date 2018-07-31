# -*- coding: utf-8 -*-
"""
Created on Fri May 18 11:21:03 2018

Plots guassian function used in simulations. 

@author: fionnlagh
"""
import numpy as np
import matplotlib.pyplot as plt

fontsize = 18
linewidth = 5
markersize = 200

unit_time = 8.5873143947605932  # s
cm_2_km = 1e-5
A = 5e6*cm_2_km  # km
Unit_l = 1e8
jet_w = 350.0/2.0/cm_2_km/Unit_l  # no units
jet_plot = [-jet_w,jet_w]

# x information
domain_nx1 = 200  # number of grid pts
AMR_lvl = 3
domain_nx1_max = domain_nx1*(2**AMR_lvl)  # max res
xprobmin1 = -4.0
xprobmax1 = 4.0
x0 = (xprobmax1-abs(xprobmin1))/2
deltax = jet_w/3.0  # @ Contians jet in guass dist as 3 sigma guass dist is zero.

# create grid pts
x = np.linspace(xprobmin1, xprobmax1, domain_nx1_max)

# y information
domain_nx2 = 20
xprobmin2 = 0.0
xprobmax2 = 6.0
y0 = 0.0  # (xprobmax2-xprobmin2)/2
deltay = (xprobmax2-xprobmin2)/domain_nx2
# create grid points
y = np.linspace(xprobmin2, xprobmax2, domain_nx2)

dist_test_x = A*np.exp(-((x-x0)/deltax)**2)
dist_test_y = A*np.exp(-((y-y0)/deltay)**2)

plt.xlim(-0.2,0.2)
plt.rc('xtick', labelsize=fontsize)
plt.rc('ytick', labelsize=fontsize)
plt.xlabel('x [Mm]', fontsize=fontsize)
plt.ylabel('Velocity [km s-1]', fontsize=fontsize)
plt.scatter(jet_plot, np.zeros(len(jet_plot)),alpha=0.9, s=markersize, zorder=1, color='r')
plt.plot(x, dist_test_x, linewidth=linewidth, zorder=-1)
plt.grid()
plt.show()
