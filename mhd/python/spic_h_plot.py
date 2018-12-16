# -*- coding: utf-8 -*-
"""
Created on Fri May 18 11:21:03 2018

Plots guassian function used in simulations.

@author: fionnlagh
"""
import numpy as np
import matplotlib.pyplot as plt
import math

V = np.arange(3, 7)*10
B = np.arange(3, 9)*10
indexs = np.meshgrid(V,B)
# data = np.zeros((len(V), len(B)))

r0 = [3.38E+08, 3.60E+08, 3.94E+08, 4.16E+08, 4.39E+08, 4.50E+08]
r1 = [4.61E+08, 4.95E+08, 5.29E+08, 5.63E+08, 5.96E+08, 6.08E+08]
r2 = [5.42E+08, 6.75E+08, 7.20E+08, 7.54E+08, 7.88E+08, 7.99E+08]
r3 = [8.33E+08, 8.55E+08, 9.11E+08, 9.79E+08, 1.00E+09, 1.00E+09]
data = r0
data = np.vstack((data, r1))
data = np.vstack((data, r2))
data = np.vstack((data, r3))

fontsize = 18
plt.rc('xtick', labelsize=fontsize)
plt.rc('ytick', labelsize=fontsize)
plt.ylabel('Height [Mm]', fontsize=fontsize)

T = True
F = False
V_vs_data = F
offset = 0.1
dV = V[1]-V[0]
dB = B[1]-B[0]
if V_vs_data is True:
    plt.xlabel('Velocity [km s-1]', fontsize=fontsize)
    plt.xlim(V[0]-offset*dV, V[-1]+offset*dV)
    plot_array = np.zeros((len(B), len(V)))
    plt.ylim(3, 10.5)
    for i in range(len(B)):
        plt.plot(V, data[:, i]*1e-8, '-o', linewidth=2, label=str(B[i])+'G')
    plt.legend(loc=2)
else:
    plt.xlabel('Magnetic Field Strength [G]', fontsize=fontsize)
    plt.xlim(B[0]-offset*dB, B[-1]+offset*dB)    
    plot_array = np.zeros((len(V), len(B)))
    plt.ylim(3, 11.5)
    for i in range(len(V)):
        plt.plot(B, data[i]*1e-8, '-o', linewidth=2, label=str(V[i])+'km s-1')
    plt.legend(loc=2)

plt.grid()
plt.show()

plt.show()
