import numpy as np
import glob
import img2vid as i2v
import os

mini_path = '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet'
sav_loc = mini_path+'/vids'

B = ['30', '40', '50', '60', '70', '80']
V = ['30', '40', '50', '60']
fps = 7

for ii in range(len(B)):
    for jj in range(len(V)):
        prefix = 'jet_B'+str(B[ii])+'_V'+str(V[jj])
        file_sav_name = 'jet_B'+str(B[ii])+'_V'+str(V[jj])+'_.'
        i2v.image2video(filepath=sav_loc+'/jet_B'+str(B[ii])+'_V'+str(V[jj]), prefix=file_sav_name, in_extension='png', 
                        output_name=prefix+'video2', out_extension='avi', 
                        fps=fps, n_loops=1, delete_images=True, 
                        delete_old_videos=True, res=1080, overlay=False, cover_page=False)
