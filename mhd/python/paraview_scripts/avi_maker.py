
import img2vid as i2v
save_directory = '/home/fionnlagh/work/AMR_code/mhd/python/paraview_scripts/testing/images'
prefix = 'test.'
fps = 4
res = (1920, 1080)

i2v.image2video(filepath=save_directory, prefix=prefix, in_extension='png', 
                output_name=prefix+'video', out_extension='avi', 
                fps=fps, n_loops=1, delete_images=True, 
                delete_old_videos=True, res=1080, overlay=False, cover_page=False)