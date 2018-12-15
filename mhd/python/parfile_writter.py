# -*- coding: utf-8 -*-
"""
Created on Mon Aug 13 16:36:22 2018

@author: fionnlagh
"""
import os
from collections import OrderedDict

path = "/home/smp16fm/work/AMR_code/solar_atmos/solar_atmosphere_2.5D/sharc_parfiles"
os.chdir(path)

# create parfile info
filelist = OrderedDict([('base_filename=', "'test'"),
            ('saveprim=', '.T.'),
            ('autoconvert=', '.T.'),
            ('convert_type=', "'vtuBCCmpi'"),
            ('slice_type=', "'csv'"),
            ('nwauxio=', '6')])

savelist = OrderedDict([('itsave(1,1)=', '0'),
            ('itsave(1,2)=', '0'),
            ('dtsave_log=', '0.25d0'),
            ('dtsave_dat=', '0.25d0'),
            ('dtsave(3)=', '0.25d0 !<= slice output saved'),
            ('nslices=', '1 !<number of slices'),
            ('slicedir(1)=', '1    !<= data is taken perp to dim select'),
            ('slicecoord(1)=', '0.0d0 !<= pt of origin of slice.')])
stoplist = OrderedDict([('dtmin=', '1.D-7',
            'time_max=', '50.0d0')])
methodlist = OrderedDict([('time_integrator=', "'twostep'"),
              ('flux_scheme=', '20*'+"'tvdlf'"),
              ('typepred1=', '20*'+"'hancock'"),
              ('typelimited=', "'predictor'"),
              ('small_pressure=', '1.0d-8'),
              ('small_density=', '1.0d-14')])
boundlist = OrderedDict([('nghostcells=', '2',
             ('typeboundary_min1 =', "'symm','asymm','symm','asymm','symm','asymm','symm','asymm'"),
             ('typeboundary_max1 =', "'symm','asymm','symm','asymm','symm','asymm','symm','asymm'"),
             ('typeboundary_min2 =', '8*'+"'special'"),
             ('typeboundary_max2 =', '8*'+"'special'")])
meshlist = OrderedDict([('refine_max_level=', '4'),
            ('refine_criterion=', '3'),
            ('refine_threshold=', '20*0.2d0'),
            ('derefine_ratio=', '20*0.15d0'),
            ('w_refine_weight(1)=', '0.4d0'),
            ('w_refine_weight(6)=', '0.2d0'),
            ('w_refine_weight(7)=', '0.2d0'),
            ('w_refine_weight(8)=', '0.2d0'),
            ('block_nx1=', '10'),
            ('block_nx2=', '10'),
            ('domain_nx1=', '300'),
            ('domain_nx2=', '100'),
            ('iprob=', '1'),
            ('xprobmin1=', '-4.d0'),
            ('xprobmax1=', '4.d0'),
            ('xprobmin2=', '0.d0'),
            ('xprobmax2=', '90.d0')])
paramlist = OrderedDict([('typecourant=', "'maxsum'"),
             'courantpar=', '0.8d0')])
mhd_list = OrderedDict([('mhd_thermal_conduction=', '.false.'),
            ('mhd_radiative_cooling=', '.false.'),
            ('mhd_gravity=', '.true.'),
            ('typedivbfix=', "'linde'"),
            ('Busr=', '20.d0'),
            ('B0field=', '.false.'),
            ('boundary_divbfix(3)=', '.false.'),
            ('mhd_n_tracer=', '1')])
# Custom entries
atmos_list = OrderedDict([('Te_profile=', "'C7'"),
              ('npts=', '8000')])
my_parameters = OrderedDict([('amp =', '50.0d0 !km s-1'),
                 ('B_strength =', '60.0d0 !G'),
                 ('jet_time =', '240.0d0 !s'),
                 ('alpha_val =', '0.0d0'),
                 ('tilt_pc =', '0.0d0')])
my_switches = OrderedDict([('tanh_profile=', '.T.'),
               ('c7_profile =', '.F.'),
               ('integrate =', '.T.'),
               ('derivative =', '.F.'),
               ('derivative_2 =', '.F.'),
               ('driver=', '.F.'),
               ('driver_random =', '.F.'),
               ('driver_kuz =', '.F.'),
               ('driver_injetion =', '.F.'),
               ('jet_cont =', '.F.'),
               ('jet_switch_on_off =', '.T.'),
               ('jet_skewed_guass =', '.F.')])
template = OrderedDict([('&filelist', filelist),
            ('&savelist', savelist),
            ('&stoplist', stoplist),
            ('&methodlist', methodlist),
            ('&boundlist', boundlist),
            ('&meshlist', meshlist),
            ('&paramlist', paramlist),
            ('&mhd_list', mhd_list),
            ('&atmos_list', atmos_list),
            ('&my_parameters', my_parameters),
            ('&my_switches', my_switches)])


def parfile_creation(master_dir, par_path, par_name, sav_path, sav_loc, template):
    # Creates dirs for parfile and qsubs
    if os.path.isdir(master_dir+par_path) is False:
        os.makedirs(master_dir+par_path)
    # creates dirs for the location of the save sim op
    if os.path.isdir(sav_loc+sav_path) is False:
        os.makedirs(sav_loc+sav_path)
    template['&filelist']['base_filename='] = "'"+sav_loc+sav_path+'/'+par_name+"_'"
    file = open(master_dir+par_path+'/'+par_name+'.par', 'w')
    for i in template:
        file.write(str(i)+"\n")
        for j in template[i]:
            file.write(str(j)+str(template[i][j])+"\n")
        file.write('/ \n \n')
    file.close()


def submitter_creation(master_dir, par_path, par_name, nb_cores, email, rmem, run_time):
    file_sub = open(master_dir+par_path+'/sub_'+par_name, 'w')
    file_sub.write('#!/bin/bash \n#$ -l h_rt=' + run_time+'\n' +
                   '#$ -pe mpi ' + nb_cores + '\n#$ -l rmem='
                   + rmem +'\n#$ -m bea' + '\n#$ -M ' + email + '\n#$ -j y'
                   + '\n\nmodule load mpi/openmpi/2.1.1/gcc-6.2' 
                   +'\n\nmpirun -np ' +nb_cores+' '+master_dir+'/amrvac -i .'+par_path+'/'+par_name+'.par')
    file_sub.close() 


def bash_writter(qsub_names_list):
    file_bash = open('multi_qsub', 'w')
    file_bash.write('#!/bin/bash\n')
    for i_list in range(len(qsub_names_list)):
        file_bash.write('qsub '+qsub_names_list[i_list]+' &\n')
    file_bash.close()
# NOTE: need to make sure there is a / infront of shared
sav_loc = 'shared/mhd_jet1/User/smp16fm/sims/jet2'
# corosponding submitters.
run_time = '40:00:00'
nb_cores = '24'
rmem = '2G'
email = 'fmackenziedover1@sheffield.ac.uk'

B = ['30', '40', '50', '60', '70', '80']
V = ['30', '40', '50', '60']

qsub_names_list = []
for bi in range(len(B)):
    template['&my_parameters']['B_strength ='] = B[bi]+'d0'
    for vi in range(len(V)):
        par_name = 'jetB'+str(B[bi])+'V'+str(V[vi])
        par_path =  '/fresh/B'+str(B[bi])'/V'+str(V[vi])
        sav_path = /fresh/B'+str(B[bi])'/V'+str(V[vi])
        template['&my_parameters']['amp ='] = V[vi]+'d0'        
        parfile_creation(master_dir, par_path, par_name, sav_path, sav_loc, template)
        submitter_creation(master_dir, par_path, par_name, nb_cores, email, rmem, run_time)
        qsub_names_list.append('.'+par_path+'/sub_'+par_name)

bash_writter(qsub_names_list)
