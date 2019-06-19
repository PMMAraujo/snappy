import os
import argparse
import subprocess

FOLDERS = ['trees','aligned','blast']
FILES = ['subtype_results.csv', 'report_subtype_results.csv',
         'all_aligned.fasta']

def clean_folders(path):
    """Folder cleaner.

    This function iterates over the folders defined in 'FOLDERS' list deleting
    the files inside them and after that the folder itself.

	Args:
        path (str): String representing the path to a folder.

	Returns:
        This function does not return.
	"""
    if os.path.exists(path):
        in_folder = os.listdir(path)
        for entry in in_folder:
            os.remove(f'{path}/{entry}')
        os.rmdir(path)

def clean_files(path):
    """File cleaner.

    This function iterates over the files defined in 'FILES' list deleteing
    them.

	Args:
        path (str): String representing the path to a file.

	Returns:
        This function does not return.
	"""
    if os.path.exists(path):
        os.remove(path)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-w","--what", required = True, help="What to delete.")
    args = parser.parse_args()
    if args.what == 'all':
        for folder in FOLDERS:
            clean_folders(folder)
        for file in FILES:
            clean_files(file)

    elif args.what == 'inter':
        for folder in FOLDERS:
            clean_folders(folder)

    elif args.what == 'c_all':
        subprocess.call('tar -zcf all.tar.gz *', shell=True)
        for folder in FOLDERS:
            clean_folders(folder)
        for file in FILES:
            clean_files(file)

    elif args.what == 'c_inter':
        if os.path.exists('trees'):
            subprocess.call("tar -zcf intermidiate_files.tar.gz split \
            aligned trees", shell=True)
            for folder in FOLDERS:
                clean_folders(folder)
        else:
            subprocess.call("tar -zcf intermidiate_files.tar.gz split aligned",
            shell=True)
            for folder in FOLDERS:
                clean_folders(folder)

    else:
        print('ERROR: Invalid parameter passed to the "cleanner.py" script.')
