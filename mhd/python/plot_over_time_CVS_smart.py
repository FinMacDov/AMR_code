import csv
import matplotlib.pyplot as plt
from matplotlib import interactive
from matplotlib import rcParams, cycler
from scipy.interpolate import interp1d
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import img2vid as i2v
import glob
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

linewidth = 4

colour_map = F 
line_plot_evo = F
vs_cs = F

plot_Te_v2 = T
movie = T

mass_density_sum = F

plasma_beta = F

sound_speed = F

img_save = T

location = '/media/fionnlagh/W7_backup/c7test/longrun'
save_loc = '/home/fionnlagh/work/AMR_code/mhd/python/image_testing/'

step = 2
fps = 4
fontsize = 18
labelsize = 14
B_cgs = 50  # G
Guass_2_telsa = 1e-4
B_SI = B_cgs*Guass_2_telsa  # tesla
miu0_si = 1.2566370614e-6  # H m^-1
# cgs
runit = 2.3416704877999998E-015
punit = 0.31754922400000002
tunit = 1e6
gamma = 5/3
vunit = 11645084.295622544
timeunit = 155.658*2#2.14683

# location = '/home/fionnlagh/work/dat/mhd/2d/c7_relax_run_csv'
#location = '/home/fionnlagh/sim_data/jet_simple/slice1'


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
tick_clip = []
pt = []
rhot = []
yt = []
Tet = []

spt = 0
ept = spt+40 #len(fname)
NBpts = 800

tick = np.arange(spt, ept)*timeunit
for sim_time in range(spt, ept, step):
    tick_clip.append(sim_time*timeunit)
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
Te_ulimit = np.floor(np.max(Tet)+0.5*np.max(Tet))
Te_llimit = np.floor(np.min(Tet)-0.5*np.min(Tet))
mach = np.sqrt(gamma*(np.array(pt)/np.array(rhot)))
v2t_cs_array = v2t/mach

if colour_map:
    H, Tunit = np.meshgrid(yt, tick_clip)
    # https://matplotlib.org/examples/color/colormaps_reference.html
    fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, figsize=(18, 20),
                                        dpi=80, facecolor='w', edgecolor='k')
    if vs_cs == T:
        v_cs_limit = np.floor(np.max(np.abs(v2t_cs_array)))
        levels1 = MaxNLocator(nbins=500).tick_values(-v_cs_limit, v_cs_limit)
        cmap1 = plt.get_cmap('seismic')
        norm1 = BoundaryNorm(levels1, ncolors=cmap1.N)
        cf1 = ax1.contourf(Tunit, H, v2t_cs_array, cmap=cmap1, levels=levels1)
        fig.colorbar(cf1, ax=ax1, label='$\frac{v_y}{c_s}$ [$\mathrm{m\;s^{-1}}$]')  # v2 units
    else:
        levels1 = MaxNLocator(nbins=500).tick_values(-vlimit, vlimit)
        cmap1 = plt.get_cmap('seismic')
        norm1 = BoundaryNorm(levels1, ncolors=cmap1.N)
        cf1 = ax1.contourf(Tunit, H, v2t, cmap=cmap1, levels=levels1)
        fig.colorbar(cf1, ax=ax1, label='$v_y$ [$\mathrm{m\;s^{-1}}$]')  # v2 units

    cmap2 = plt.get_cmap('hot')
    levels2 = MaxNLocator(nbins=100).tick_values(np.min(rhot), np.max(rhot))
    norm2 = BoundaryNorm(levels2, ncolors=cmap2.N)

    cf2 = ax2.contourf(Tunit, H, rhot, cmap=cmap2, levels=levels2)
    fig.colorbar(cf2, ax=ax2, label=r'$\rho$ [$\mathrm{kg\;m^{-3}}$]')
    ax2.set_ylabel('$y$ [$\mathrm{Mm}$]', fontsize=14)

    cmap3 = plt.get_cmap('coolwarm')
    levels3 = MaxNLocator(nbins=100).tick_values(np.min(Tet), np.max(Tet))
    norm3 = BoundaryNorm(levels3, ncolors=cmap3.N)

    cf3 = ax3.contourf(Tunit, H, Tet, cmap=cmap3, levels=levels3)
    fig.colorbar(cf3, ax=ax3, label='Te [$\mathrm{K}$]')
    ax3.set_xlabel('Time [s]', fontsize=14)
    if img_save:
        fig.savefig(save_loc+'colour_map_vy_rho_Te')
    plt.close(fig)
if line_plot_evo:
    spacer = 0
    zero_pts = np.zeros(len(v2t))
    # Selects every nth row
    LIM = [0, 10]
    v2t_cs = v2t_cs_array[LIM[0]::LIM[-1]].T
    # Make color gradient
    N = len(v2t_cs.T)
    cmap = plt.cm.coolwarm
    rcParams['axes.prop_cycle'] = cycler(color=cmap(np.linspace(0, 1, N)))

    SPACER = 5
    gaps = np.linspace(0, len(v2t_cs.T)*SPACER, len(v2t_cs.T))
    line = np.full(len(gaps), 7)
    zero_line = np.full((len(yt), len(gaps)), gaps)

    fig2, ax2 = plt.subplots(figsize=(18, 20))
    fontsize = 16

    plt.plot(yt, v2t_cs+gaps)
    plt.plot(yt, zero_line, linestyle='--')

    labels = np.chararray(N, unicode=True, itemsize=5)
    yt_label_pos = yt[LIM[0]::LIM[-1]].T
    labels[:] = ' '
    labels[::4] = (str(-SPACER))
    labels[1::4] = ('0.0')
    labels[2::4] = (str(SPACER))
    plt.yticks(gaps, labels, fontsize=fontsize)
    plt.rc('xtick', labelsize=fontsize)
    time_stamps = tick[LIM[0]::LIM[-1]]
    text_pos_1 = np.full(len(v2t_cs.T), 9)
    text_pos_2 = zero_line[0]+1
    for i in range(len(text_pos_1)):
        ax2.text(text_pos_1[i], text_pos_2[i],
                 str(round(time_stamps[i], -1))+' s', fontsize=fontsize)
    plt.xlabel('y [Mm]', fontsize=fontsize)
    plt.ylabel('v2/cs', fontsize=fontsize)
    plt.xlim(0, 10)
if plot_Te_v2:
    fig, (ax1, ax2) = plt.subplots(2, 1)
    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    for i in range(len(Tet)):
        ax1.set_ylabel('Te [$\mathrm{K}$]')
        ax2.set_ylabel('$v_y$ [$\mathrm{m\;s^{-1}}$]')
        ax1.semilogy(yt, Tet[i], color='r', linewidth=linewidth)
        ax1.semilogy(yt, Tet[0], color='grey', linewidth=linewidth,
                     linestyle='dashed')
        ax1.set_ylim(Te_llimit, Te_ulimit)

        ax2.plot(yt, v2t[i], linewidth=linewidth)
        ax2.set_ylim(-vlimit-vlimit*0.1, vlimit+vlimit*0.1)
        textstr = '%.3f' % round(tick_clip[i]/60**2, 3)+' hrs'
        ax2.text(0.75, 0.1, textstr, transform=ax2.transAxes, 
                 fontsize=fontsize,bbox=props)
        if img_save:
            fig.savefig(save_loc+'/evo/'+'Te_vy_evo'+"{0:04}".format(i),  dpi=180)
        ax1.clear()
        ax2.clear()
    plt.close(fig)
    if movie:
        master_dir = save_loc+'evo'
        name = 'Te_vy_evo'
        total_files = glob.glob(master_dir+'/'+name+'*.png')
        
        hsize = 586
        wsize = 782
        sav_path = '/shared/mhd_jet1/User/smp16fm/sims/atmos/c7/vids'
        
        total_files = sorted(total_files)
        
        prefix = name
        file_sav_name = name
        sav_loc = master_dir
        
        i2v.image2video(filepath=sav_loc, prefix=file_sav_name, in_extension='png', 
                        output_name=file_sav_name+'video', out_extension='avi', 
                        fps=fps, n_loops=1, delete_images=True, 
                        delete_old_videos=True, res=1080, overlay=False, cover_page=False)

if mass_density_sum:
    total_mass_den = np.asarray(rhot).sum(axis=1)
    plt.plot(tick_clip/60**2,total_mass_den, linewidth=linewidth)
    plt.ylabel('Total mass density [$\mathrm{kg\;m^{-3}}$]',  fontsize=fontsize)
    plt.xlabel('Time [hrs]',  fontsize=fontsize)
    plt.tick_params(labelsize=labelsize)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    if img_save:
        plt.savefig(save_loc+'total_mass_density')
    plt.close()
if plasma_beta:
    plt.semilogy(yt, pt[0]/(B_SI**2/2*miu0_si), linewidth=linewidth)
    plt.ylabel(r'$\beta$',  fontsize=fontsize)
    plt.xlabel('$y$ [$\mathrm{Mm}$]',  fontsize=fontsize)
    plt.tick_params(labelsize=labelsize)
    if img_save:
        plt.savefig(save_loc+'plasma_beta')
    plt.close()
if sound_speed:
    plt.semilogy(yt, mach[0], linewidth=linewidth)
    plt.ylim(mach[0].min()-0.1*mach[0].min(), mach[0,-1]+0.1*mach[0,-1])
    plt.ylabel('$c_S$ [$\mathrm{m\;s^{-1}}$]',  fontsize=fontsize)
    plt.xlabel('$y$ [$\mathrm{Mm}$]',  fontsize=fontsize)
    plt.tick_params(labelsize=labelsize)
    if img_save:
        plt.savefig(save_loc+'sound_speed')
    plt.close()
