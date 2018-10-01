# -*- coding: utf-8 -*-
"""
Created on Mon Aug 13 16:36:22 2018

@author: fionnlagh
python3 /home/smp16fm/work/AMR_code/mhd/python/parfile_writter_P_scan.py 
"""
import os
path = "/home/smp16fm/work/AMR_code/solar_atmos/solar_atmosphere_2.5D/sharc_parfiles"
os.chdir(path)

# create parfile info
filelist = {'base_filename=': "'test'",
            'saveprim=': '.T.',
            'autoconvert=': '.T.',
            'convert_type=': "'vtuBCCmpi'",
            'slice_type=': "'csv'",
            'nwauxio=': '6'}

savelist = {'itsave(1,1)=': '0',
            'itsave(1,2)=': '0',
            'dtsave_log=': '0.25d0',
            'dtsave_dat=': '0.25d0',
            'dtsave(3)=': '0.25d0 !<= slice output saved',
            'nslices=': '1 !<number of slices',
            'slicedir(1)=': '1    !<= data is taken perp to dim select',
            'slicecoord(1)=': '0.0d0 !<= pt of origin of slice.'}
stoplist = {'dtmin=': '1.D-7',
            'time_max=': '50.0d0'}
methodlist = {'time_integrator=': "'twostep'",
              'flux_scheme=': '20*'+"'tvdlf'",
              'typepred1=': '20*'+"'hancock'",
              'typelimited=': "'predictor'",
              'small_pressure=': '1.0d-8',
              'small_density=': '1.0d-14'}
boundlist = {'nghostcells=': '2',
             'typeboundary_min1 =': "'symm','asymm','symm','asymm','symm','asymm','symm','asymm'",
             'typeboundary_max1 =': "'symm','asymm','symm','asymm','symm','asymm','symm','asymm'",
             'typeboundary_min2 =': '8*'+"'special'",
             'typeboundary_max2 =': '8*'+"'special'"}
meshlist = {'refine_max_level=': '4',
            'refine_criterion=': '3',
            'refine_threshold=': '20*0.2d0',
            'derefine_ratio=': '20*0.15d0',
            'w_refine_weight(1)=': '0.4d0',
            'w_refine_weight(6)=': '0.2d0',
            'w_refine_weight(7)=': '0.2d0',
            'w_refine_weight(8)=': '0.2d0',
            'block_nx1=': '10',
            'block_nx2=': '10',
            'domain_nx1=': '300',
            'domain_nx2=': '100',
            'iprob=': '1',
            'xprobmin1=': '-4.d0',
            'xprobmax1=': '4.d0',
            'xprobmin2=': '0.d0',
            'xprobmax2=': '90.d0'}
paramlist = {'typecourant=': "'maxsum'",
             'courantpar=': '0.8d0'}
mhd_list = {'mhd_thermal_conduction=': '.false.',
            'mhd_radiative_cooling=': '.false.',
            'mhd_gravity=': '.true.',
            'typedivbfix=': "'linde'",
            'Busr=': '20.d0',
            'B0field=': '.false.',
            'boundary_divbfix(3)=': '.false.',
            'mhd_n_tracer=': '1'}
# Custom entries
atmos_list = {'Te_profile=': "'C7'",
              'npts=': '8000'}
my_parameters = {'amp =': '50.0d0 !km s-1',
                 'B_strength =': '60.0d0 !G',
                 'jet_time =': '300.0d0 !s',
                 'alpha_val =': '0.0d0',
                 'tilt_pc =': '0.0d0'}
my_switches = {'tanh_profile=': '.T.',
               'c7_profile =': '.F.',
               'integrate =': '.T.',
               'derivative =': '.F.',
               'derivative_2 =': '.F.',
               'driver=': '.F.',
               'driver_random =': '.F.',
               'driver_kuz =': '.F.',
               'driver_injetion =': '.F.',
               'jet_cont =': '.T.',
               'jet_switch_on_off =': '.F.',
               'jet_skewed_guass =': '.F.'}
template = {'&filelist': filelist,
            '&savelist': savelist,
            '&stoplist': stoplist,
            '&methodlist': methodlist,
            '&boundlist': boundlist,
            '&meshlist': meshlist,
            '&paramlist': paramlist,
            '&mhd_list': mhd_list,
            '&atmos_list': atmos_list,
            '&my_parameters': my_parameters,
            '&my_switches': my_switches}

# corosponding submitters.
run_time = '10:00:00'
nb_cores = '24'
rmem = '2G'
email = 'fmackenziedover1@sheffield.ac.uk'

# Input of changing vaibles.
B = ['50']
V = ['50']
P = ['160','200','240','280','320','360']
save_path = '/shared/mhd_jet1/User/smp16fm/sims/jet/'
qsub_run_path = 'sharc_parfiles/'
excute_path = '/home/smp16fm/work/AMR_code/solar_atmos/solar_atmosphere_2.5D/amrvac'

file_bash = open('multi_qsub', 'w')
file_bash.write('#!/bin/bash\n')
for bi in range(len(B)):
    for vi in range(len(V)):
        for pi in range(len(P)):
            if os.path.isdir('B'+B[bi]+'/V'+V[vi]+'/P'+P[pi]+'/') is False:
                os.makedirs('B'+B[bi]+'/V'+V[vi]+'/P'+P[pi]+'/')
    
            if os.path.isdir(save_path+'/jet_B'+B[bi]+'_V'+V[vi]+'P'+P[pi]) is False:
                os.makedirs(save_path+'/jet_B'+B[bi]+'_V'+V[vi]+'P'+P[pi])
    
            file = open('B'+B[bi]+'/V'+V[vi]+'/P'+P[pi]+'/B'+B[bi]+'_V'+V[vi]+'_P'+P[pi]+'.par', 'w')
            # create submitters
            file_sub = open('B'+B[bi]+'/V'+V[vi]+'/P'+P[pi]+'/sub_B'+B[bi]+'_V'+V[vi]+'_P'+P[pi], 'w')
            file_sub.write('#!/bin/bash \n#$ -l h_rt=' + run_time+'\n' +
                           '#$ -pe mpi ' + nb_cores + '\n#$ -l rmem='
                           + rmem +'\n#$ -m bea' + '\n#$ -M ' + email + '\n#$ -j y'
                           + '\n\nmodule load mpi/openmpi/2.1.1/gcc-6.2' 
                           +'\n\nmpirun -np ' +nb_cores+' '+excute_path+' -i '+qsub_run_path+'B'+B[bi]+'/V'+V[vi]+'/P'+P[pi]+'/B'+B[bi]+'_V'+V[vi]+'_P'+P[pi]+'.par')
            file_bash.write('qsub sharc_parfiles/B'+B[bi]+'/V'+V[vi]+'/P'+P[pi]+'/sub_B'+B[bi]+'_V'+V[vi]+'_P'+P[pi]+' &\n')
            #creates file path
            template['&filelist']['base_filename='] = "'"+save_path+'jet_B'+B[bi]+'_V'+V[vi]+'P'+P[pi]+'/jet_B'+B[bi]+'_V'+V[vi]+'_P'+P[pi]+"_'"
            template['&my_parameters']['B_strength ='] = B[bi]+'.0d0  !G'
            template['&my_parameters']['amp ='] = V[vi]+'.0d0  !km s-1'
            template['&my_parameters']['jet_time ='] = P[pi]+'.0d0 !s'
            for i in template:
                file.write(str(i)+"\n")
                for j in template[i]:
                    file.write(str(j)+str(template[i][j])+"\n")
                file.write('/ \n \n')
            file.close()
file_bash.close()
