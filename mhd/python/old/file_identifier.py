def userinput(read_csv = False):
    """
    Finds the total number of FITS files in a given directory.

    Example call:

    >> [in]: from dtp import userinput
    >> [in]: userinput()
    >> [out]: Where are your files stored?
    >> [in]: Downloads
    >> [out]: resulting filepath
    etc..

    Note: Your files must be stored in the same directory as this script.

    """

    import os

    def listdir_fullpath(d):
            return [os.path.join(d, f) for f in os.listdir(d)]


    if read_csv == True:
        location = input("\nWhere are your files stored? ")
        global file_path
        file_path = os.path.abspath(location)
        global total_files
        total_files = listdir_fullpath(location)
        total_files = sorted(total_files)
        if os.path.exists(file_path):
            print("\nDirectory exists.. ")
            # checking
            counter = 0
            for file in total_files:
                if file.endswith("s"):
                    counter += 1
        if len(total_files) == counter:
            print("\nTotal number of FITS files in this directory is %s. \n" %len(total_files)) 
        else:
            raise Exception("\nError! Not all files in the directory are FITS.")
    else:
        raise Exception("\nPlease set read_fits/read_slits to True/False.")