import img2vid as i2v

master_dir = '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/v_jet/hd/2D'
fps = 7

jet_angle = ['0.1','0.5','1.0','5.0','10.0','15.0','20.0','25.0','30.0']
for ii in range(len(jet_angle)):
        sav_loc = master_dir+'/vids/jet_'+jet_angle[ii]
        file_sav_name = 'jet_a'+jet_angle[ii]+'.'
        i2v.image2video(filepath=sav_loc, prefix=file_sav_name, in_extension='png', 
                        output_name=file_sav_name+'video', out_extension='avi', 
                        fps=fps, n_loops=1, delete_images=True, 
                        delete_old_videos=True, res=1080, overlay=False, cover_page=False)
