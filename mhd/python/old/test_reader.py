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
import os


def listdir_fullpath(d):
    return [os.path.join(d, f) for f in os.listdir(d)]

location = '/home/fionnlagh/work/dat/mhd/2d/c7_relax_run_csv'
global file_path
file_path = os.path.abspath(location)
global total_files
total_files = listdir_fullpath(location)
total_files = sorted(total_files)

fnames = []
for file in total_files:
    if file.endswith("csv"):
        fnames.append(file)
