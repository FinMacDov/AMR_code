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
F = False
T = True

VALC = F
C7 = F
McWhirter = F
VALMC = T
Testing = T


if VALC is True:
    os.chdir('/home/fionnlagh/work/AMR_code/mhd/python/VAL3C_data')
    df = pd.read_csv('formatted_valc3c.csv')
    df = df.sort_values(by=['z_h'])
    df2 = pd.read_csv('valc_intpr.csv')

# endif

if C7 is True:
    os.chdir('/home/fionnlagh/work/AMR_code/mhd/python/c7_data')
    df = pd.read_csv('C7_raw_csv.csv')
    df = df.sort_values(by=['z_h'])
    df2 = pd.read_csv('c7_intpr.csv')
# endif

if McWhirter is True:
    os.chdir('/home/fionnlagh/work/AMR_code/mhd/python/McWhirter_data')
    df = pd.read_csv('formatted_McWhirter.csv')
    df = df.sort_values(by=['z_h'])
    df2 = pd.read_csv('McWhirter_intpr.csv')

if VALMC is True:
    os.chdir('/home/fionnlagh/work/AMR_code/mhd/python')
    dfVALC = pd.read_csv('VAL3C_data/formatted_valc3c.csv')
    dfMc = pd.read_csv('McWhirter_data/formatted_McWhirter.csv')
    os.chdir('/home/fionnlagh/work/AMR_code/mhd/python/VALMC_data')
    dfVALC = dfVALC.sort_values(by=['z_h'])
    dfMc = dfMc.sort_values(by=['z_h'])
    df2 = pd.read_csv('VALMC_intpr.csv')
    df_alt = pd.read_csv('VALMC_intpr_edited.csv')
    z_alt = df_alt.z_h.values
    Te_alt = df_alt.Te.values
    plt.semilogy(z_alt, Te_alt, linewidth='3', color='gold')



if VALMC is False:
    varname = list(df)
    if VALC is True:
        rho_orig = df.rho.values
    # endif
    Te_orig = df.Te.values
    Te_intrp = df2.Te.values

    if McWhirter is True:
        z_orig = df.z_h.values
        z_intrp = df2.z_h.values
    else:
        z_orig = df.z_h.values*km_2_cm
        z_intrp = df2.z_h.values
     # endif
# endif
else:
    Te_orig_VALC = dfVALC.Te.values
    Te_orig_Mc = dfMc.Te.values

    z_orig_VALC = dfVALC.z_h.values*km_2_cm
    z_orig_Mc = dfMc.z_h.values

    Te_orig = []
    z_orig = []
    for i in range(len(Te_orig_VALC)):
        Te_orig.append(Te_orig_VALC[i])
        z_orig.append(z_orig_VALC[i])
    for j in range(len(Te_orig_Mc)):
        Te_orig.append(Te_orig_Mc[j])
        z_orig.append(z_orig_Mc[j])
    z_intrp = df2.z_h.values
    Te_intrp = df2.Te.values

plt.xlabel('Height [cm]')
plt.ylabel('Te [K]')


#plt.plot(z_orig, Te_orig, linewidth='3')
#plt.plot(z_intrp, Te_intrp, linestyle='--', linewidth='3')    
#plt.semilogy(z_orig, Te_orig,linestyle='--', linewidth='3',zorder=2)
plt.xlim(0, 1e9)
plt.ylim(1e3, 1e7)
#plt.semilogy(z_intrp, Te_intrp, linewidth='3',zorder=1)
#plt.scatter(z_intrp, Te_intrp, color='red', zorder=1)
#plt.scatter(z_orig, Te_orig, color='black', zorder=2)
if Testing is True:
    p1 = 0
    p2 = -15 #-11
    p11 = 1
    p22 = -1
#    plt.xlim(2.109e8, 2.110e8)
#    plt.ylim(1.3e4, 1.4e4)
    plt.plot(z_orig_VALC[p1:p2], Te_orig_VALC[p1:p2], linewidth='3', color='red')
    plt.plot(z_orig_Mc[p11:p22]+2.2975e7, Te_orig_Mc[p11:p22], linewidth='3', color='red')

plt.show()
