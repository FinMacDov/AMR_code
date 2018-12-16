import img2vid as i2v
import glob
import os
from PIL import Image

def resizeImage(hsize,wsize,sav_path,name):
    img = Image.open(sav_path+'/'+name)
    img = img.resize((wsize,hsize), Image.ANTIALIAS)
    return img.save(sav_path+'/'+name) 
    
master_dir = '/home/fionnlagh/work/AMR_code/mhd/python/image_testing'
fps = 80
name = 'fb'

total_files = glob.glob(master_dir+'/'+name+'*.png')

hsize = 586
wsize = 782
sav_path = '/shared/mhd_jet1/User/smp16fm/sims/atmos/c7/vids'
## Renaming files to allow for ordering
#for i in range(len(total_files)):
#    os.rename(master_dir+'/'+name+str(i)+'.png', master_dir+'/'+name+"{0:04}".format(i)+'.png')
#    resizeImage(hsize,wsize,sav_path,'/test'+"{0:04}".format(i)+'.png')

total_files = sorted(total_files)

prefix = name
file_sav_name = name
sav_loc = master_dir

i2v.image2video(filepath=sav_loc, prefix=file_sav_name, in_extension='png', 
                output_name=file_sav_name+'video', out_extension='avi', 
                fps=fps, n_loops=1, delete_images=True, 
                delete_old_videos=True, res=1080, overlay=False, cover_page=False)
