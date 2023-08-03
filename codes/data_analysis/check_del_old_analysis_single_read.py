
import os;
import sys;
sys.path.insert(0, '/home/madhu/work/codes')

import string;
import glob;
import time
import copy
from multiprocessing.pool import ThreadPool, Pool
from multiprocessing import cpu_count

from tqdm import tqdm
from main import check_endWith, check_create_dir, check_conda_env

# import datetime
from datetime import datetime

import h5py
import numpy as np
import multiprocessing

from collections import defaultdict

import tempfile
import subprocess

import re;

'''
Madhu:
I have created this python script to check if all the old analysis have been properly deleted from the SINGLE read FAST5 files. 
'''


'''
(base) qian@lambda-dual-qgenlab:~/Projects/NanoporeRNA/epinano$ h5ls -r Curlcake_non/1204670/1/GXB01170_20180726_FAH85615_GA10000_mux_scan_RNAAB089716_23886_read_12_ch_193_strand.fast5
/                        Group
/Analyses                Group
/Analyses/Basecall_1D_000 Group
/Analyses/Basecall_1D_000/BaseCalled_template Group
/Analyses/Basecall_1D_000/BaseCalled_template/Events Dataset {1803}
/Analyses/Basecall_1D_000/BaseCalled_template/Fastq Dataset {SCALAR}
/Analyses/Basecall_1D_000/Configuration Group
/Analyses/Basecall_1D_000/Configuration/basecall_1d Group
/Analyses/Basecall_1D_000/Summary Group
/Analyses/Basecall_1D_000/Summary/basecall_1d_template Group
/Analyses/Calibration_Strand_Detection_000 Group
/Analyses/Calibration_Strand_Detection_000/Configuration Group
/Analyses/Calibration_Strand_Detection_000/Configuration/calib_detector Group
/Analyses/Calibration_Strand_Detection_000/Summary Group
/Analyses/Calibration_Strand_Detection_000/Summary/calibration_strand_template Group
/Analyses/Segmentation_000 Group
/Analyses/Segmentation_000/Configuration Group
/Analyses/Segmentation_000/Configuration/stall_removal Group
/Analyses/Segmentation_000/Summary Group
/Analyses/Segmentation_000/Summary/segmentation Group
/PreviousReadInfo        Group
/Raw                     Group
/Raw/Reads               Group
/Raw/Reads/Read_12       Group
/Raw/Reads/Read_12/Signal Dataset {27050/Inf}
/UniqueGlobalKey         Group
/UniqueGlobalKey/channel_id Group
/UniqueGlobalKey/context_tags Group
/UniqueGlobalKey/tracking_id Group


*** The FAST5 files with multiple read also look like: 

(ACONDA) madhu@qgenlabgpu:~/datasets/download/p1_RNA_mod_RNA_seq_xPore/extracted_data$ h5ls -r HEK293T-KD-Ctrl-rep1/fast5/fast5_pass/FAN49781_pass_b26f7f0c_268.fast5 | head
/                        Group
/read_0002d9a7-61cf-4785-84d9-71b3d8528ce4 Group
/read_0002d9a7-61cf-4785-84d9-71b3d8528ce4/Analyses Group
/read_0002d9a7-61cf-4785-84d9-71b3d8528ce4/Analyses/Basecall_1D_000 Group
/read_0002d9a7-61cf-4785-84d9-71b3d8528ce4/Analyses/Basecall_1D_000/BaseCalled_template Group
/read_0002d9a7-61cf-4785-84d9-71b3d8528ce4/Analyses/Basecall_1D_000/BaseCalled_template/Fastq Dataset {SCALAR}
/read_0002d9a7-61cf-4785-84d9-71b3d8528ce4/Analyses/Basecall_1D_000/Summary Group
/read_0002d9a7-61cf-4785-84d9-71b3d8528ce4/Analyses/Basecall_1D_000/Summary/basecall_1d_template Group
/read_0002d9a7-61cf-4785-84d9-71b3d8528ce4/Analyses/Segmentation_000 Group
/read_0002d9a7-61cf-4785-84d9-71b3d8528ce4/Analyses/Segmentation_000/Summary Group

'''

def do_check(args):
   # t0 = datetime.now()
   a_fast5_file = args[0]
   found = False
   try: 
      f = h5py.File(a_fast5_file, 'r')
      # print(list(f.keys())[0:5])
      for a_key in list(f.keys()):
         if 'Analyses' in a_key:
            print('This file contains OLD Analysis: ', a_fast5_file)
            found = True
         elif 'read_' in a_key:
            print('This file has multiple reads: ', a_fast5_file)
            found = True
         # else:
         #    print('')
            # print('This file is fine: ', a_fast5_file)

   except Exception as e:
        print('a_fast5_file: ', a_fast5_file)
        print('Exception:', e)
   
   sys.stdout.flush()
   return(found, a_fast5_file)



if __name__=='__main__':
   if len(sys.argv)<2:
      print("Usage: {} {}".format(sys.argv[0], "basefolder", "sub_folder_level"))
      print("Example:")
      print("       python {} {} {} {}".format(sys.argv[0], "/home/madhu/datasets/download/m6A_RNA_modification_native_RNA_seq/GSM3528751/extracted_data"))
#       print("       python {} {} {} {}".format(sys.argv[0], "Curlcake_non"))
      sys.exit(1)

   startTime = datetime.now()
   current_time = startTime.strftime("%Y/%m/%d at %H:%M:%S")
   print('\n\nThe script starting at: ' + str(current_time), ' \n\n' )
   sys.stdout.flush()
                                                                                                                                                                        
   basefolder = sys.argv[1]
   # basefolder = '/home/madhu/datasets/download/p5_decoding_epitranscriptional_native_RNA_seq/single_read_extracted_data/IVT_m7G_fast5'
   basefolder = check_endWith(basefolder)
   
   
   all_fast5_files = [os.path.join(root, name)
             for root, dirs, files in os.walk(basefolder)
             for name in files
             if name.endswith((".fast5"))]
   
   found_old_analysis = False
   
   print('Time now is: ', str(datetime.now()))
   print('It took: ', str(datetime.now() - startTime), ' to get ', len(all_fast5_files) , ' fast5 files inside: ', basefolder)
   sys.stdout.flush()


   # a_fast5_file = '/home/madhu/datasets/download/p5_decoding_epitranscriptional_native_RNA_seq/extracted_data/HS_SAEC_fast5/HS_SAEC_fast5/HS_SAEC_319.fast5'
   # f = h5py.File(a_fast5_file, 'r')
   # print(list(f.keys())[0:5])
   # # if 'Analyses' in list(f.keys()) or 'read_' in list(f.keys()):
   # for a_key in list(f.keys()):
   #    if 'Analyses' in list(f.keys()) or 'read_' in list(f.keys()):
   #       print('This file contains OLD Analysis or is a multiple read: ', a_fast5_file)
   #       found_old_analysis = True
   # else:
   #    print('This file is fine: ', a_fast5_file)




   # print('There are: ', len(all_fast5_files), ' fast5 files inside: ', basefolder)
   
   inputs = zip(all_fast5_files)

   ## Use the maximum CPU possible 
   cpus = cpu_count()
   n_jobs = min(cpus, len(all_fast5_files))

   # n_jobs = 20

   print('n_jobs: ', n_jobs)
   ## I am using multiple Processors as it is faster than threads. 
   results = tqdm(Pool(n_jobs).imap_unordered(do_check, inputs), total=len(all_fast5_files))
   # results = ThreadPool(n_jobs).imap_unordered(do_check, inputs)
   for result in results:
      print()
      # if result[0] == True: 
      #    print("\n\nXXXXXXX There are files with OLD Analysis or multiple reads. XXXXXXX\n\n")
      #    print('For: ', result[1])
      # else:
      #    print('\n\nAll Good.\n\n')

   print("\nDone!")
   executionTime = (datetime.now() - startTime)
   current_time = datetime.now().strftime("%Y/%m/%d at %H:%M:%S")
   print('\n\nThe script completed at: ' + str(current_time))
   print('Execution time: ' + str(executionTime), ' \n\n')




