# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 10:53:24 2018

@author: fin
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import csv
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import numpy as np
import os
        
y   = []
rho = []
p   = []
Te  = []

os.chdir('/home/fionnlagh/work/AMR_code')
with open('atmos_data/c7/csv/C7.csv') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
         y.append(row['Y']) 
         rho.append(row['rho'])
         p.append(row['p'])
         Te.append(row['Te'])

rho = np.array(rho)
rho = rho.astype(np.float)
p = np.array(p)
p = p.astype(np.float)
Te = np.array(Te)
Te = Te.astype(np.float)
y = np.array(y)
y = y.astype(np.float) 

pnew = interp1d(y, p, kind='linear') 
rhonew = interp1d(y, rho, kind='linear')#SI
Tenew = interp1d(y, Te, kind='linear') #SI
ynew = np.linspace(y[0],y[-1],num=8000, endpoint=True)
Te_interp = Tenew(ynew)
rho_interp = rhonew(ynew)
p_interp = pnew(ynew)

grav = -274.0 # m s-1
solrad = 6.95e8 # m

ggrid = grav*(solrad/(solrad+ynew))**2

gradp = np.gradient(p_interp)
rho_grad = gradp/ggrid 

fb_grad = -gradp + rho_interp*ggrid 

pa = np.empty(len(Te_interp))
ra = np.empty(len(Te_interp))
ra[0] = rho_interp[0] 
pa[0] = ra[0]*Te_interp[0] 
inT = ggrid[0]/Te_interp[0]
for j in range(1,len(Te_interp)):
    inT=inT+(ggrid[j]/Te_interp[j]+ggrid[j-1]/Te_interp[j-1])*0.5
    pa[j]=pa[0]*np.exp(inT*ynew[1])
    ra[j]=pa[j]/Te_interp[j]

gradpa = np.gradient(pa, edge_order=2)
fba = gradpa-ra*ggrid 

maxnb = 3
ra_list = np.zeros((maxnb, len(Te_interp)))
ra_list[0] = ra
pa_list = np.zeros((maxnb, len(Te_interp)))
pa_list[0] = pa
fb_list = np.zeros((maxnb, len(Te_interp)))
fb_list[0] = fba
gradpa_list = np.zeros((maxnb, len(Te_interp)))
gradpa_list[0] = gradpa
for j in range(1,maxnb):
      ra_list[j] = gradpa_list[j-1]/ggrid + fb_list[j-1]
      pa_list[j] = ra_list[j]*Te_interp 
      gradpa_list[j] = np.gradient(pa_list[j], edge_order=2)
      fb_list[j] = gradpa_list[j]-ra*ggrid
      
F = False
T = True

plot_te = F
plot_te_sim = F
plot_fb = F
plot_p = T

#----------------------------------------------
# Multiplot
#----------------------------------------------
if plot_te is True:
    fig, ax1 = plt.subplots()
    ax1.plot(y, rho, 'b-')
    ax1.plot(ynew, rho_grad, 'g-')
    ax1.set_xlabel('Y')
    # Make the y-axis label, ticks and tick labels match the line color.
    ax1.set_ylabel('rho', color='b')
    ax1.tick_params('y', colors='black')
    plt.yscale('log')

    ax2 = ax1.twinx()
    ax2.plot(y, Te, 'r-')
    ax2.plot(ynew,Te_interp, 'g-')
    ax2.set_ylabel('Te', color='r')
    ax2.tick_params('y', colors='black')
    plt.yscale('log')

    fig.tight_layout()
    plt.show()        
#----------------------------------------------
#end        
#----------------------------------------------
#----------------------------------------------
# Multiplot
#----------------------------------------------
if plot_te_sim is True:
    fig, ax1 = plt.subplots()
    ax1.plot(y, rho, 'b-')
    ax1.plot(ynew, rho_grad, 'g-')
    ax1.set_xlabel('Y')
    # Make the y-axis label, ticks and tick labels match the line color.
    ax1.set_ylabel('rho', color='b')
    ax1.tick_params('y', colors='black')
    plt.yscale('log')

    ax2 = ax1.twinx()
    ax2.plot(y, Te, 'r-')
    ax2.plot(ynew,Te_interp, 'g-')
    ax2.set_ylabel('Te', color='r')
    ax2.tick_params('y', colors='black')
    plt.yscale('log')

    fig.tight_layout()
    plt.show()        
#----------------------------------------------
#end        
#----------------------------------------------
if plot_fb is True:
    fig, ax3 = plt.subplots()
    ax3.plot(ynew, fba, 'b-')
    for j in range(0,maxnb):
        ax3.plot(ynew, fb_list[j])
    ax3.set_xlabel('Y')
    # Make the y-axis label, ticks and tick labels match the line color.
    ax3.set_ylabel('fb', color='b')
    ax3.tick_params('y', colors='black')
#    ax3.set_ylim([np.min(fba), np.max(fba)])
    
    fig.tight_layout()
    plt.show()        
#----------------------------------------------
# Multiplot
#----------------------------------------------
if plot_p is True:
    fig, ax3 = plt.subplots()
    ax3.plot(ynew, pa, 'b-')
    for j in range(0,maxnb):
        ax3.plot(ynew, pa_list[j])
    ax3.set_xlabel('Y')
    # Make the y-axis label, ticks and tick labels match the line color.
    ax3.set_ylabel('fb', color='b')
    ax3.tick_params('y', colors='black')
#    ax3.set_ylim([np.min(fba), np.max(fba)])
    
    fig.tight_layout()
    plt.show()   

