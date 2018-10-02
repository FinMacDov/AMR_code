# -*- coding: utf-8 -*-
"""
Created on Mon Aug 13 16:36:22 2018

@author: fionnlagh
"""
import os
path = "/home/smp16fm/work/AMR_code/solar_atmos/solar_atmosphere_2.5D/sharc_parfiles"
os.chdir(path)

# create parfile info
filelist = {'base_filename=': "'test'",
            'saveprim=': '.T.',
            'autoconvert=': '.T.',
            'convert_type=': "'vtuBmpi'",
            'slice_type=': "'csv'",
            'nwauxio=': '6'}

savelist = {'itsave(1,1)=': '0',
            'itsave(1,2)=': '0',
            'ditsave_log=': '100',
            'dtsave_dat=': '1.d-2',
            'dtsave(3)=': '0.25d0 !<= slice output saved',
            'slicedir(1)=': '1    !<= data is taken perp to dim select',
            'slicecoord(1)=': '0.0d0 !<= pt of origin of slice.'}
stoplist = {'dtmin=': '1.0d-9',
            'time_max=': '0.6d0'}
methodlist = {'time_integrator=': "'threestep'",
              'flux_scheme=': '20*'+"'hllc'",
              'limiter=':"'cada3'"}
boundlist = {'nghostcells=': '2',
             'typeboundary_min1 =': '4*'+"'cont'",
             'typeboundary_max1 =': '4*'+"'cont'",
             'typeboundary_min2 =': '4*'+"'special'",
             'typeboundary_max2 =': '4*'+"'cont'"}
meshlist = {'refine_max_level=': '4',
            'refine_threshold=': '20*0.5d0',
            'derefine_ratio=': '20*1.d0/2.d0',
            'block_nx1=': '20',
            'block_nx2=': '20',
            'domain_nx1=': '240',
            'domain_nx2=': '60',
            'iprob=': '1',
            'xprobmin1=': '-0.5d0',
            'xprobmax1=': '3.d0',
            'xprobmin2=': '0.d0',
            'xprobmax2=': '1.5d0'}
paramlist = {'slowsteps=': "'1000'",
             'courantpar=': '0.6d0'}
my_parameters = {'jet_speed =': '5.0d0,
                 'plasma_beta =': '100.0d0 ',
                 'tilt_deg =': '0.0d0',
                 'jet_width=': '0.05d0',
                 'driv_tanh=': '.T.',
                 'driv_gaussian':'.F.',
                 'orig_setup': '.F.'}
template = {'&filelist': filelist,
            '&savelist': savelist,
            '&stoplist': stoplist,
            '&methodlist': methodlist,
            '&boundlist': boundlist,
            '&meshlist': meshlist,
            '&paramlist': paramlist,
            '&my_parameters': my_parameters}


def parfile_creation(master_dir, par_path, par_name, template):
    if os.getcwd() is not master_dir:
        os.chdir(master_dir)
    if os.path.isdir(master_dir+par_path) is False:
        os.makedirs(par_path)
    template['&filelist']['base_filename='] = "'"+par_path+'/'+par_name+"_'"
    file = open(master_dir+par_path+'/'+par_name+'.par', 'w')
    for i in template:
        file.write(str(i)+"\n")
        for j in template[i]:
            file.write(str(j)+str(template[i][j])+"\n")
        file.write('/ \n \n')
    file.close()


def submitter_creation(master_dir, par_path, par_name, nb_cores, email, rmem, run_time):
    file_sub = open(par_path+'/sub_'+par_name, 'w')
    file_sub.write('#!/bin/bash \n#$ -l h_rt=' + run_time+'\n' +
                   '#$ -pe mpi ' + nb_cores + '\n#$ -l rmem='
                   + rmem +'\n#$ -m bea' + '\n#$ -M ' + email + '\n#$ -j y'
                   + '\n\nmodule load mpi/openmpi/2.1.1/gcc-6.2' 
                   +'\n\nmpirun -np ' +nb_cores+' '+master_dir+'/amrvac -i '+par_path+par_name+'.par')
    file_sub.close()
    qsub_names = par_path+'/sub_'+par_name 


def bash_writter(qsub_names_list):
    file_bash = open('multi_qsub', 'w')
    file_bash.write('#!/bin/bash\n')
    for i_list in range(len(qsub_names_list)):
        file_bash.write('qsub '+qsub_names_list[i_list]+' &\n')
    file_bash.close()


master_dir = '/home/smp16fm/work/AMR_code/v_jets/hd_jets'
# corosponding submitters.
run_time = '10:00:00'
nb_cores = '24'
rmem = '2G'
email = 'fmackenziedover1@sheffield.ac.uk'

jet_angle = ['0.0']#,'0.1','0.5','1','5','10','15','20','25','30']

qsub_names_list = []
for jj in range(len(jet_angle)):
    par_name = 'HD_jet_M'+template['&my_parameters']['jet_speed ='][0:1]+'_a'+jet_angle
    par_path =  '/parfiles/jet_a'+jet_angle[jj]
    template['&my_parameters']['tilt_deg ='] = jet_angle[jj]+'.0d0'
    parfile_creation(master_dir, par_path, sav_loc, par_name, template)
    submitter_creation(master_dir, par_path, par_name, nb_cores, email, rmem, run_time)
    qsub_names_list.append(qsub_names)

bash_writter(qsub_names_list)


