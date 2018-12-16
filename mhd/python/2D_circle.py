# -*- coding: utf-8 -*-
"""
Created on Fri May 18 11:21:03 2018

Plots guassian function used in simulations.

@author: fionnlagh
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def calc_e(p, hd_gamma, rho, v_x, v_y):
    e = (p/(hd_gamma-1.0))+0.5*rho*((v_x)**2+(v_y)**2)
    return e
    
def calc_e_mag(p, hd_gamma, rho, v_x, v_y,b1,b2):
    e_mag = (p/(hd_gamma-1.0))+0.5*rho*(v_x**2+v_y**2)+0.5*(b1**2+b2**2)
    return e_mag

T = True
F = False

jet_speed = F
jet_rho = F
jet_e = T
jet_e_mag = F
jet_m = F
line_plot = F

fontsize = 18
linewidth = 5
markersize = 200

refine_max_level = 4
domain_nx1 = 240
domain_nx2 = 60
xprobmin1 = -0.5
xprobmax1 = 3.0
xprobmin2 = 0.0
xprobmax2 = 1.5
res_y = 1.5/domain_nx2
nghost = 2.0

domain_nx1_max = domain_nx1*(2**refine_max_level)  # max res
domain_nx2_max = domain_nx2*(2**refine_max_level)  # max res

hd_gamma = 1.4
rhoj = hd_gamma
vj = 5.0
eta = 3.0
rho_b = rhoj/eta
p = 1.0
tilt_deg = 0.0
tilt = (np.pi/180.0)*tilt_deg 
beta = 0.5
b1 = 0.0 
b2 = np.sqrt(2.0/beta)

x = np.linspace(xprobmin1, xprobmax1, domain_nx1_max)  # spatial pos x
y = np.linspace(nghost*res_y-xprobmin2, xprobmax2, domain_nx2_max)  # spatial pos y

X,Y = np.meshgrid(x, y) # spatial grid

jet_w = 0.5
# centre pos of circle
xc = 0.0
yc = 0.0
# converts grid into polar co-orndinates
rcloud = (X-xc)**2+(Y-yc)**2
#condition to give tanh bc
trwidth = jet_w  # transition width
PoI = 0.5*(-jet_w+(-jet_w+trwidth))  # point of inflection
smoothness = trwidth*0.1
v1_tanh_con = np.sin(tilt)*vj*0.5*(1.0-np.tanh((np.sqrt(rcloud)+PoI)/smoothness))
v2_tanh_con = np.cos(tilt)*vj*0.5*(1.0-np.tanh((np.sqrt(rcloud)+PoI)/smoothness))
rho_tanh_con = rho_b+(rhoj-rho_b)*0.5*(1.0-np.tanh((np.sqrt(rcloud)+PoI)/smoothness))
e = calc_e(p, hd_gamma, rho_tanh_con, v1_tanh_con, v2_tanh_con)
e_mag = calc_e_mag(p, hd_gamma, rho_tanh_con, v1_tanh_con, v2_tanh_con,b1,b2) 
v_tot = np.sqrt(v1_tanh_con**2+v2_tanh_con**2)

if jet_speed:
    fig, ax1 = plt.subplots()
    p1 = ax1.pcolormesh(X, Y, v_tanh_con, cmap=cm.plasma, vmin=v_tanh_con.min(), vmax=v_tanh_con.max())
    plt.xlim(x[0], x[-1])
    plt.ylim(y[0], y[-1])
    cb = fig.colorbar(p1)
if jet_rho:
    fig, ax2 = plt.subplots()
    p2 = ax2.pcolormesh(X, Y, rho_tanh_con, cmap=cm.plasma, vmin=rho_tanh_con.min(), vmax=rho_tanh_con.max()) 
    plt.xlim(x[0], x[-1])
    plt.ylim(y[0], y[-1])
    cb = fig.colorbar(p2)
if jet_e:
    fig, ax3 = plt.subplots()
    p3 = ax3.pcolormesh(X, Y, e, cmap=cm.plasma, vmin=e.min(), vmax=e.max())   
    plt.xlim(x[0], x[-1])
    plt.ylim(y[0], y[-1])
    cb = fig.colorbar(p3)
if jet_m:
    fig, ax4 = plt.subplots()
    p4 = ax4.pcolormesh(X, Y, rho_tanh_con*v1_tanh_con, cmap=cm.plasma, vmin=((rho_tanh_con*v1_tanh_con)
).min(), vmax=(rho_tanh_con*v1_tanh_con).max())   
    plt.xlim(x[0], x[-1])
    plt.ylim(y[0], y[-1])
    cb = fig.colorbar(p4)
    fig, ax5 = plt.subplots()
    p5 = ax5.pcolormesh(X, Y, rho_tanh_con*v2_tanh_con, cmap=cm.plasma, vmin=((rho_tanh_con*v2_tanh_con)
).min(), vmax=(rho_tanh_con*v2_tanh_con).max())   
    plt.xlim(x[0], x[-1])
    plt.ylim(y[0], y[-1])
    cb = fig.colorbar(p5)
if jet_e_mag:
    fig, ax6 = plt.subplots()
    p6 = ax6.pcolormesh(X, Y, e_mag, cmap=cm.plasma, vmin=e_mag.min(), vmax=e_mag.max())   
    plt.xlim(x[0], x[-1])
    plt.ylim(y[0], y[-1])
    cb = fig.colorbar(p6)
if line_plot:
    slice_rho_0 = v2_tanh_con[np.where(X==x[np.argmin(abs(x))])]
    y_loc = Y[np.where(X==x[np.argmin(abs(x))])]
    plt.plot(y_loc,slice_rho_0, linewidth = linewidth, markersize = markersize)

#fig, ax = plt.subplots()
#
#im = plt.imshow(rho_grid, cmap=cm.plasma, vmin=rho_grid.min(), vmax=rho_grid.max(), extent=[x[0], x[-1], y[0], y[-1]])
#im.set_interpolation('bilinear')
#
#cb = fig.colorbar(im)


