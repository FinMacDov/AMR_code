# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 10:53:24 2018

@author: fin
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os


cm_2_km = 1e-5
km_2_cm = 1e5
Mm_to_cm = 1e8
F = False
T = True

# Raw data
os.chdir('/home/fionnlagh/work/AMR_code/mhd/python/c7_data')
df = pd.read_csv('C7_raw_csv.csv')
df = df.sort_values(by=['z_h'])

Te_orig = df.Te.values
z_orig = df.z_h.values*km_2_cm

# Amrvac data
os.chdir('/home/fionnlagh/work/dat/mhd/2d/jet_test')
df_amr = pd.read_csv('c7_test_2__d1_x+0.00D+00_n0000.csv')

tunit = 1e6
punit = 0.27612976
runit = 1.6726217770000001e-15

p_amr = df_amr.p.values
rho_amr = df_amr.rho.values
Te_amr = ((p_amr/punit)/(rho_amr/runit))*tunit
z_amr = df_amr.Y.values*Mm_to_cm

plt.xlabel('Height [cm]')
plt.ylabel('Te [K]')

#plt.plot(z_orig, Te_orig, linewidth='3')   
plt.semilogy(z_orig, Te_orig, linestyle='--', linewidth='3', zorder=2)
plt.semilogy(z_amr, Te_amr, linewidth='3', zorder=1)
#plt.plot(z_amr, Te_amr, linewidth='3', zorder=1)
#plt.scatter(z_orig, Te_orig, color='black', zorder=2)
plt.show()
