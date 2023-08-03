import os


import sys
sys.path.insert(0, '/home/madhu/datasets/codes')

from main import check_endWith, check_create_dir, check_conda_env


'''
This script will check all the FAST5 ffiles and directly put them under the given folder. 
'''


base_folder = '/home/madhu/datasets/download/p3_m6A_RNA_modification_native_RNA_seq/extracted_data/WT1_RNAAA023484'

FAST5_files = [os.path.join(root, name)
             for root, dirs, files in os.walk(base_folder)
             for name in files
             if name.endswith((".fast5"))]

print(FAST5_files[1:5])

### I have decided to move files manually, as it is more safe .... 