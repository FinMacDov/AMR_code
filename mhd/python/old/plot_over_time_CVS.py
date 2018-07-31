# -*- coding: utf-8 -*-
"""
Created on Wed May  9 13:12:34 2018

@author: fionnlagh
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import csv
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import numpy as np
import os

os.chdir('/home/fionnlagh/work/dat/mhd/2d/relax_run')
# set file name
fname = 'relax_tanh_big_box_d1_x+0.00D+00_n'
timeunit = 171.746
# cgs
runit = 2.3416704877999998E-015
punit = 0.31754922400000002
tunit = 1e6
gamma = 5/3
vunit = 11645084.295622544
v2t = []
yt = []
tick = []
pt = []
rhot = []
yt = []
Tet = []


for time in range(0, 5):
    print(time)
    if time in range(0, 10):
        nbname = str(0)+str(0)+str(0)+str(time)
    if time in range(11, 100):
        nbname = str(0)+str(0)+str(time)
    if time in range(101, 1000):
        nbname = str(0)+str(time)
    with open(fname+nbname+'.csv') as csvfile:
        reader = csv.DictReader(csvfile)
        y = []
        rho = []
        p = []
        v2 = []
        for row in reader:
            rho.append(row['rho'])
            p.append(row['p'])
            v2.append(row['v2'])
            y.append(row['Y'])
        rho = np.array(rho)
        rho = rho.astype(np.float)
        p = np.array(p)
        p = p.astype(np.float)
        Te = np.divide(p/punit, rho/runit)*tunit
        y = np.array(y)
        y = y.astype(np.float)
        v2 = np.array(v2)
        v2 = v2.astype(np.float)
        cs = np.sqrt(gamma*vunit*0.01*Te/tunit)  # SI
        pnew = interp1d(y, p, kind='cubic')
        rhonew = interp1d(y, rho*1000, kind='cubic')  # SI
#        rhonew = interp1d(y, rho/runit, kind='cubic')  # no dim
        v2new = interp1d(y, v2*0.01, kind='cubic')  # SI
        Tenew = interp1d(y, Te, kind='cubic')  # SI
        ynew = np.linspace(y[0], y[-1], num=500, endpoint=True)
        v2_interp = v2new(ynew)
        rho_interp = rhonew(ynew)
        p_interp = pnew(ynew)
        Te_interp = Tenew(ynew)
        tick.append(time*timeunit)
        yt.append(ynew)
        v2t.append(v2_interp)
        pt.append(p_interp)
        rhot.append(rho_interp)
        Tet.append(Te_interp)

vlimit = np.floor(np.max(np.abs(v2t)))  # v2 units
# vlimit = round(np.max(np.abs(v2t)),3)  # unitless
Tunit, H = np.meshgrid(ynew, tick)

# https://matplotlib.org/examples/color/colormaps_reference.html

levels1 = MaxNLocator(nbins=500).tick_values(-vlimit, vlimit)
cmap1 = plt.get_cmap('seismic')
norm1 = BoundaryNorm(levels1, ncolors=cmap1.N)

fig, (ax1, ax2, ax3) = plt.subplots(nrows=3)

cf1 = ax1.contourf(H, Tunit, v2t, cmap=cmap1, levels=levels1)
fig.colorbar(cf1, ax=ax1, label='v2 [m s-1]')  # v2 units
# fig.colorbar(cf, ax=ax1,label='v2/cs')# unitless
# ax1.set_ylabel('y [Mm]', fontsize=14)
# ax1.set_xlabel('Time [s]', fontsize=14)

cmap2 = plt.get_cmap('hot')
levels2 = MaxNLocator(nbins=100).tick_values(np.min(rhot), np.max(rhot))
norm2 = BoundaryNorm(levels2, ncolors=cmap2.N)

cf2 = ax2.contourf(H, Tunit, rhot, cmap=cmap2, levels=levels2)
fig.colorbar(cf2, ax=ax2, label='rho [kg m-3]')
ax2.set_ylabel('y [Mm]', fontsize=14)
#ax2.set_xlabel('Time [s]', fontsize=14)    

#cmap2 = plt.get_cmap('hot')
#levels2 = MaxNLocator(nbins=100).tick_values(np.min(pt),np.max(pt))
#norm2 = BoundaryNorm(levels2, ncolors=cmap2.N)
#
#cf2 = ax2.contourf(H, Tunit, pt, cmap=cmap2,levels=levels2)
#fig.colorbar(cf2, ax=ax2, label='p [dyne cm-2]')
#ax2.set_ylabel('y [Mm]', fontsize=14)
##ax2.set_xlabel('Time [s]', fontsize=14)    


cmap3 = plt.get_cmap('coolwarm')
levels3 = MaxNLocator(nbins=100).tick_values(np.min(Tet), np.max(Tet))
#levels3 = MaxNLocator(nbins=100).tick_values(4.5e3,8e5)
norm3 = BoundaryNorm(levels3, ncolors=cmap3.N)

cf3 = ax3.contourf(H, Tunit, Tet, cmap=cmap3, levels=levels3)
fig.colorbar(cf3, ax=ax3, label='Te [K]')
#ax3.set_ylabel('y [Mm]', fontsize=14)
ax3.set_xlabel('Time [s]', fontsize=14)
plt.show()
#plt.colorbar(label='T [K]')#, ticks=np.linspace(-150,150,10, endpoint=True));            