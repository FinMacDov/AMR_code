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
import pandas as pd
import numpy as np
from decimal import Decimal
import os

# Switches
VALC = True
C7 = False

# Loads correct data
if VALC is True:
    os.chdir('/home/fionnlagh/work/AMR_code/mhd/python/VAL3C_data')
    df = pd.read_csv('formatted_valc3c.csv')
    df = df.sort_values(by=['z_h'])
# endif

if C7 is True:
    os.chdir('/home/fionnlagh/work/AMR_code/mhd/python/c7_data')
    df = pd.read_csv('C7_raw_csv.csv')
    df = df[::-1]
# endif

varname = list(df)

km_to_cm = 1e5

test = list(df.Te.values)

test_from = []
for k in range(len(test)):
    test[k] = ("{:.5E}".format(test[k]))

remains = int(np.floor((len(test)/6+len(test))/6))

for l in range(1, remains+1):
    test.insert(l*6-1, '&')

test[0] = str.replace(str(test[0]), 'e', 'd')
test[0] = str.replace(str(test[0]), 'E', 'd')
test[0] = str.replace(str(test[0]), '+', '')
test[0] = str.replace(str(test[0]), 'd0', 'd')

file = open('testfile.txt', 'w')
file.write(str(test[0])+', ')
# formats data for FORTRAN table
for j in range(1, len(test)):
    markerj = (j+1) % 6
    test[j] = str.replace(str(test[j]), 'e', 'd')
    test[j] = str.replace(str(test[j]), 'E', 'd')
    test[j] = str.replace(str(test[j]), '+', '')
    test[j] = str.replace(str(test[j]), 'd0', 'd')
    if j == (len(test)-1):
        file.write(str(test[j]))
    elif markerj == 0:
        file.write(str(test[j])+' \n                          ')
    else:
        file.write(str(test[j])+', ')
file.close()
