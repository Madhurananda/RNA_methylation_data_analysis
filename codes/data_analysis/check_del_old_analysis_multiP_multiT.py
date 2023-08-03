
import os;
import sys;
import string;
import glob;
import time
import copy

from tqdm import tqdm
from multiprocessing.pool import Pool, ThreadPool
import threading

import concurrent.futures


import h5py
import numpy as np
import multiprocessing

from collections import defaultdict

import tempfile
import subprocess

import re;

'''
Madhu:
I have created this python script to check if all the old analysis have been properly deleted. 
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

'''

def check_analyses(all_fast5_files):
   for a_fast5_file in tqdm(all_fast5_files):
      f = h5py.File(a_fast5_file, 'r')
      if 'Analyses' in list(f.keys()) or '/read_' in f:
         print('This file contains OLD Analysis or is a multiple read: ', a_fast5_file)
         found_old_analysis = True
      else:
         print('This file is fine: ', a_fast5_file)


if __name__=='__main__':
   if len(sys.argv)<2:
      print("Usage: {} {}".format(sys.argv[0], "basefolder", "sub_folder_level"))
      print("Example:")
      print("       python {} {} {} {}".format(sys.argv[0], "/home/madhu/datasets/download/m6A_RNA_modification_native_RNA_seq/GSM3528751/extracted_data"))
#       print("       python {} {} {} {}".format(sys.argv[0], "Curlcake_non"))
      sys.exit(1)
   
   basefolder = sys.argv[1]
   
   if basefolder.endswith('/'):
      print('')
   else:
      basefolder+='/'
   
   
   all_fast5_files = [os.path.join(root, name)
             for root, dirs, files in os.walk(basefolder)
             for name in files
             if name.endswith((".fast5"))]
   
   ## Trying automatic processes
   with concurrent.futures.ProcessPoolExecutor() as executor:
      tqdm([executor.submit(check_analyses, all_fast5_files) for _ in range(20)])
      
#    ## Trying automatic threads (working fine)
#    with concurrent.futures.ThreadPoolExecutor() as executor:
#       [executor.submit(check_analyses, all_fast5_files) for _ in range(20)]
   
#    ## Trying manual Threads (working fine)
#    threads = []
#    for _ in range(20):
#       t = threading.Thread(target=check_analyses, args=(all_fast5_files,))
#       t.start()
#       threads.append(t)

#    for thread in threads:
#       thread.join()

        
        
        
        
        

#    found_old_analysis = False
#    pmanager = multiprocessing.Manager();
#    handlers = []
#    h5files_Q = pmanager.Queue();
#    for h5f in all_fast5_files:
#       print('\n\nThe file put here is: ', h5f)
#       h5files_Q.put( h5f )
   
#    for hid in range( 20 ):
#       p = multiprocessing.Process(target=check_analyses, args=(h5files_Q,));
#       p.start();
#       handlers.append(p);
   
#    while any(p.is_alive() for p in handlers):
#       try:
#          pass;
#       except:
#          time.sleep(1);
#          continue;

#    print("Done!")
    
    
#    with ThreadPool(processes=20) as p:
# #       p.imap(check_analyses, all_fast5_files)
#       tqdm(p.imap(check_analyses, all_fast5_files), total = len(all_fast5_files))
# #       p.close()
# #       p.join()
    
#    with Pool(processes=20) as p: 
#       results = tqdm(p.imap_unordered(check_analyses, all_fast5_files), total = len(all_fast5_files))
# #       # shutdown the process pool
# #       p.close()
# #       # wait for all issued task to complete
# #       p.join()
   
#    print('all_fast5_files: ', len(all_fast5_files))
   
    
    
#    for a_fast5_file in tqdm(all_fast5_files):
#       f = h5py.File(a_fast5_file, 'r')
#       if 'Analyses' in list(f.keys()):
#          print('This file contains OLD Analysis: ', a_fast5_file)
#          found_old_analysis = True
#       else:
#          print('This file is fine: ', a_fast5_file)


#    if found_old_analysis == True: 
#       print("\n\nXXXXXXX There are files with OLD Analysis. XXXXXXX\n\n")
#    else:
#       print('\n\nAll Good.\n\n')

