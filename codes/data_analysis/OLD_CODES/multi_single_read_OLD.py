import os,sys
sys.path.insert(0, '/home/madhu/datasets/codes')

import glob
import tarfile
import gzip
from tqdm import tqdm
from main import check_endWith, check_create_dir, check_conda_env
import numpy as np




if __name__=='__main__':
   if len(sys.argv)<3:
      print('Usage: {} {} {}'.format('multi_single_read.py', 'basefolder', 'output_folder'))
      print("Example:")
      print('       {} {} {} {}'.format('/home/madhu/datasets/scripts/chris_scripts/untar_org.py', "/home/madhu/datasets/download/m6A_RNA_modification_native_RNA_seq/GSM3528751/original_data/extracted_data/fast5/", "/home/madhu/datasets/download/m6A_RNA_modification_native_RNA_seq/GSM3528751/splitted_extracted_data/", '2'))
      print('\n\nThe basefolders should contain the zip files. There are three types of zipped files: 1: *.fast5.gz, 2: *.fast5.tar, 3: *.tar.gz* \n\n')
#       print("       {} {} {}".format(sys.argv[0], "Curlcake_m6a", "1368282"))
#       print("       {} {} {}".format(sys.argv[0], "Curlcake_non", "1204670"))
#       print("       {} {} {}".format(sys.argv[0], "Curlcake_non", "1368288"))
      sys.exit(1)
   
   basefolder = sys.argv[1]
   output_folder = sys.argv[2]
#    zip_type = int(sys.argv[3])

   basefolder = check_endWith(basefolder)
   
   print('\n\nbasefolder:', basefolder)
   print('output_folder: ', output_folder)
   
   fast5files = [os.path.join(root, name)
             for root, dirs, files in os.walk(basefolder)
             for name in files
             if name.endswith((".fast5"))]
   
   print('The number of FAST5 files: ', len(fast5files))
   list_dir_part = []
   for a_fast5files in tqdm(fast5files):
#       print('The number of compressed zipped files: ', len(fast5files))
#       print('The entire file is : ', a_fast5files)
      a_fast5files_part = a_fast5files.split('/')[-1]
#       print('The ONLY file is: ', a_fast5files_part)
      dir_part = a_fast5files.split(a_fast5files_part)[0]
#       print('The dir part is: ', dir_part)
      list_dir_part.append(dir_part)
   
   unique_list_dir_part = np.unique(np.array(list_dir_part))
#    print('The original len: ', len(list_dir_part))
#    print('The unique lan: ', len(unique_list_dir_part))
   
   for a_dir_part in tqdm(unique_list_dir_part):
      ex_dir_part = a_dir_part.split(basefolder)[-1]
      print('The ex_dir_part part is: ', ex_dir_part)
      save_dir = output_folder + ex_dir_part
      print('The save_dir part is: ', save_dir)
      
      check_create_dir(save_dir)
      
      # if not os.path.isdir(save_dir):
      #    os.makedirs(save_dir, exist_ok=True)
      #    print('**\nThe Folder has been created: ', save_dir, '**\n')
      # else:
      #    print('The folder already exists: ', save_dir)
      
      os.system( "multi_to_single_fast5 --input {} --save_path {} --recursive".format( a_dir_part, save_dir ) )


   print('\n\n')
