import numpy as np
import glob
from pydrive.auth import GoogleAuth
from pydrive.drive import GoogleDrive
import os

mini_path = '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet'
sav_loc = mini_path+'/vids'

B = ['30', '40', '50', '60', '70', '80']
V = ['30', '40', '50', '60']
fps = 7

gauth = GoogleAuth()
gauth.LocalWebserverAuth() # Creates local webserver and auto handles authentication
#
drive = GoogleDrive(gauth)

#Create folder
folder_metadata = {'title' : 'MyFolder', 'mimeType' : 'application/vnd.google-apps.folder'}
folder = drive.CreateFile(folder_metadata)
folder.Upload()

#Get folder info and print to screen
foldertitle = folder['title']
folderid = folder['id']
print('title: %s, id: %s' % (foldertitle, folderid))

#Upload file to folder
file = drive.CreateFile({"parents": [{"kind": "drive#fileLink", "id": folderid}]})
file.SetContentFile('test.png')
file.Upload()
## Example how to create file
#file1 = drive.CreateFile({'title': 'Hello2.txt'})  # Create GoogleDriveFile instance with title 'Hello.txt'.
#file1.SetContentString('Hello World! This is for further testing') # Set content of the file from given string.
#file1.Upload()
