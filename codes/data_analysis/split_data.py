
import os;
import sys;
import string;
import glob;
import time
import copy

from tqdm import tqdm

import h5py
import numpy as np
import multiprocessing

from collections import defaultdict

import tempfile
import subprocess

import re;

'''
Madhu:
I have created this python script to split the fast5 files into 15 GB chunks ....  
'''

# calculate file size in KB, MB, GB
def convert_bytes(size):
    """ Convert bytes to KB, or MB or GB"""
    for x in ['bytes', 'KB', 'MB', 'GB', 'TB']:
        if size < 1024.0:
            return "%3.1f %s" % (size, x)
        size /= 1024.0


        
if __name__=='__main__':
   if len(sys.argv)<2:
      print("Usage: {} {}".format(sys.argv[0], "basefolder"))
      print("Example:")
      print("       python {} {} {} {}".format(sys.argv[0], "/home/madhu/datasets/download/m6A_RNA_modification_native_RNA_seq/GSM3528751/extracted_data"))
#       print("       python {} {} {} {}".format(sys.argv[0], "Curlcake_non"))
      sys.exit(1)
   
   basefolder = sys.argv[1]
   output_dir = sys.argv[2]
#    output_dir = '/home/madhu/datasets/download/m6A_RNA_modification_native_RNA_seq/GSM3528751/splitted_extracted_data'

   if basefolder.endswith('/'):
      print('')
   else:
      basefolder+='/'
   
   if output_dir.endswith('/'):
      print('')
   else:
      output_dir+='/'
   
   ## Read all the fast5 files in the basefolder 
   all_fast5_files = list(glob.iglob(basefolder+'*.fast5', recursive=False))
   
    
    
   all_fast5_files = all_fast5_files[:50000]
   cut_off_GB = 0.01
#    cut_off_GB = 15
    
    
    
    
   print('The number of files are: ', len(all_fast5_files))
   
   file_list = []
   temp_file_list = []
   
   total_size = 0
   cut_off_size = cut_off_GB*1024*1024*1024
   for a_file in tqdm(all_fast5_files):
      dir_size = os.path.getsize(a_file)
#       print('The folder: ', a_file, ' has size: ', round(dir_size/(1024*1024), 4), ' MB')
      temp_file_list.append(a_file)
      total_size += dir_size
      if total_size > cut_off_size:
         print('\n\nTotal size has been :', round((total_size/(1024*1024*1024)),4), ' GB')
         print('There are ', len(temp_file_list), ' files\n\n')
         file_list.append(temp_file_list)
         total_size = 0
         temp_file_list = []

   dir_n = 0
   # Now, move all the splitted files to seperate folders
   for i in tqdm(range(len(file_list))): 
      a_chunk_files = file_list[i]
      print('The file list of length of ', len(a_chunk_files), ' is being moved.')
      for a_c_file in a_chunk_files:
         os.system("mkdir {}".format( output_dir+'split_'+str(i+1) ));
         os.system("cp {} {}".format( a_c_file, output_dir+'split_'+str(i+1) ));
#          cp a_c_file output_dir+'split_'+str(i+1)
#    os.system("mkdir {}".format( output_dir+'split_'+str(i+1) ));
#    os.system("cp {} {}".format( a_c_file, a_c_file output_dir+'split_'+str(i+1) ));


'''
I am actually extracting each zip file (tar/gzip) for each folder. 
So, I will probably not be using this script at all. 
'''


