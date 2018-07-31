import csv
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import numpy as np
import os


def listdir_fullpath(d):
    return [os.path.join(d, f) for f in os.listdir(d)]

# cgs
runit = 2.3416704877999998E-015
punit = 0.31754922400000002
tunit = 1e6
gamma = 5/3
vunit = 11645084.295622544

rho_cgs_2_si = 1e3

location = '/home/fionnlagh/work/AMR_code/mhd/python/test'
global file_path
file_path = os.path.abspath(location)
global total_files
total_files = listdir_fullpath(location)
total_files = sorted(total_files)
timeunit = 171.746

v2t = []
yt = []
tick = []
pt = []
rhot = []
yt = []
Tet = []


with open(total_files[0]) as csvfile:
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

#plt.subplot(2, 1, 1)
#plt.semilogy(y, Te, '-', linewidth=3.0)
#plt.ylim(5.5e3, 1.5e6)
#plt.ylabel('T[k]')
#
#plt.subplot(2, 1, 2)
#plt.semilogy(y, rho*rho_cgs_2_si, '-', linewidth=3.0)
#plt.xlabel('y (Mm)')
#plt.ylabel(r'$\rho$ [kg m-3]')


fig, ax1 = plt.subplots()
fontsize = 14
linewidth = 3
color = 'red'
ax1.set_xlabel('Height (Mm)', fontsize=fontsize)
ax1.set_ylim(5.5e3, 1.5e6)
ax1.set_ylabel('Temperature [k]', fontsize=fontsize)
line1, = ax1.semilogy(y, Te, color=color, linewidth=linewidth, label='line1')
ax1.tick_params(axis='y')

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'blue'
ax2.set_ylabel(r'$\rho$ [kg m-3]', fontsize=fontsize)
line2, = ax2.semilogy(y, rho*rho_cgs_2_si, color=color, linewidth=linewidth, linestyle = '--', label = 'line2')
ax2.tick_params(axis='y')

plt.legend([line1, line2], ['Te', r'$\rho$'])

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()


plt.show()
