import csv
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import numpy as np
import os

"""
Created on Wed May  9 13:12:34 2018

@author: Fionnlagh Mackenzie Dover
"""


def listdir_fullpath(d):
    return [os.path.join(d, f) for f in os.listdir(d)]

F = False
T = True

colour_map = F
line_plot = T

# cgs
runit = 2.3416704877999998E-015
punit = 0.31754922400000002
tunit = 1e6
gamma = 5/3
vunit = 11645084.295622544
timeunit = 155.658*2

# location = '/home/fionnlagh/work/dat/mhd/2d/c7_relax_run_csv'
# location = '/media/fionnlagh/W7_backup/c7test/longrun'
location = '/media/fionnlagh/W7_backup/jet_simple/lvl4/slice1'

global file_path
file_path = os.path.abspath(location)
global total_files
total_files = listdir_fullpath(location)
total_files = sorted(total_files)
fname = []
for file in total_files:
    if file.endswith("csv"):
        fname.append(file)

v2t = []
yt = []
tick = []
pt = []
rhot = []
yt = []
Tet = []

spt = 0
ept = len(fname)
NBpts = 800

tick = np.arange(spt, ept)*timeunit
for sim_time in range(spt, ept):
    with open(fname[sim_time]) as csvfile:
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
        # This allows us to work with the data
        rho = np.array(rho).astype(np.float)
        p = np.array(p).astype(np.float)
        Te = np.divide(p/punit, rho/runit)*tunit
        y = np.array(y).astype(np.float)
        v2 = np.array(v2).astype(np.float)
        # Interpolation of the 1d data set for current time step
        p_interp = interp1d(y, p*0.1, kind='linear')  # SI
        rho_interp = interp1d(y, rho*1000, kind='linear')  # SI
        v2_interp = interp1d(y, v2*0.01, kind='linear')  # SI
        Te_interp = interp1d(y, Te, kind='linear')  # SI
        # Selects data points for interp
        yt = np.linspace(y[0], y[-1], num=NBpts, endpoint=True)
        # Collect data together
        v2t.append(v2_interp(yt))
        rhot.append(rho_interp(yt))
        pt.append(p_interp(yt))
        Tet.append(Te_interp(yt))

vlimit = np.floor(np.max(np.abs(v2t)))

if colour_map == T:
    H, Tunit = np.meshgrid(yt, tick)

    # https://matplotlib.org/examples/color/colormaps_reference.html
    levels1 = MaxNLocator(nbins=500).tick_values(-vlimit, vlimit)
    cmap1 = plt.get_cmap('seismic')
    norm1 = BoundaryNorm(levels1, ncolors=cmap1.N)

    fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, figsize=(8, 10),
                                        dpi=80, facecolor='w', edgecolor='k')

    cf1 = ax1.contourf(Tunit, H, v2t, cmap=cmap1, levels=levels1)
    fig.colorbar(cf1, ax=ax1, label='v2 [m s-1]')  # v2 units

    cmap2 = plt.get_cmap('hot')
    levels2 = MaxNLocator(nbins=100).tick_values(np.min(rhot), np.max(rhot))
    norm2 = BoundaryNorm(levels2, ncolors=cmap2.N)

    cf2 = ax2.contourf(Tunit, H, rhot, cmap=cmap2, levels=levels2)
    fig.colorbar(cf2, ax=ax2, label='rho [kg m-3]')
    ax2.set_ylabel('y [Mm]', fontsize=14)

    cmap3 = plt.get_cmap('coolwarm')
    levels3 = MaxNLocator(nbins=100).tick_values(np.min(Tet), np.max(Tet))
    norm3 = BoundaryNorm(levels3, ncolors=cmap3.N)

    cf3 = ax3.contourf(Tunit, H, Tet, cmap=cmap3, levels=levels3)
    fig.colorbar(cf3, ax=ax3, label='Te [K]')
    ax3.set_xlabel('Time [s]', fontsize=14)
    plt.show()
elif line_plot == T:
    fig1 = plt.figure(figsize=(14, 12))
    spacer = 0
    for i in range(len(v2t)):
        mach = np.sqrt(gamma*(pt[i]/rhot[i]))
        v2t_cs = v2t/mach
        v2t_cs[i] += spacer
        spacer += vlimit/np.max(mach)
        if i % 10 == 1:
            plt.plot(yt, np.full(len(v2t_cs[i]), spacer), linestyle='--')
            plt.plot(yt, v2t_cs[i])
    plt.xlabel('y [Mm]')
    plt.ylabel('v2/cs')
    plt.show()
else:
    print('I do nothing')
