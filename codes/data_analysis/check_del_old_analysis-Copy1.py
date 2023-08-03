
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

# def del_analysis( subfolderlevel, h5files_Q):
#    while not h5files_Q.empty():
#       try:
#          f5files = h5files_Q.get(block=False)
#       except:
#          break;
    
#       print('\n\nThe f5files files are: ', f5files)
#       print('The length of the f5files, ', len(f5files))
#       print('subfolderlevel is: ', subfolderlevel)  
    
#       if subfolderlevel<1:
            
# #           for f5t in f5files:
# #             print('An f5t file: ', f5t, '\n\n\n')
# #             try:
# #                 with h5py.File(f5t, 'r+') as mf5:
# #                     del mf5['/Analyses']
# #                     mf5.flush();
# #                     mf5.close();
# #             except:
# #                 mf5.close();
# #                 print("Failed to delete analysis for {}".format(f5t)) 
        
        
# #          print('The length of the f5files, ', len(f5files))
#          with h5py.File( f5files, 'r+') as mf5:
#            print('mf5 is: ', mf5)
#            for g_k in mf5.keys():
#                print('g_k is: ', g_k)
#                del mf5['/'+g_k+'/Analyses']
#            mf5.flush();
#            mf5.close();
#            print('\n\n')
        
#       else:
#          f5files = glob.glob(os.path.join( f5files, "*.fast5") )
#          print('The length of the f5files, ', len(f5files))
#          for f5t in f5files:
#             #print('An f5t file: ', f5t, '\n\n\n')
#             try:
#                 with h5py.File(f5t, 'r+') as mf5:
#                     del mf5['/Analyses']
#                     mf5.flush();
#                     mf5.close();
#             except:
#                 mf5.close();
#                 print("Failed to delete analysis for {}".format(f5t))

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
   
#    sub_level = 2
   sub_level = sys.argv[2]
   search_str = '*/'*sub_level
   ## Read all the fast5 files in the basefolder 
   all_fast5_files = list(glob.iglob(basefolder+search_str+'*.fast5', recursive=False))
   
   found_old_analysis = False
   
#    print('all_fast5_files: ', all_fast5_files)
   
   for a_fast5_file in tqdm(all_fast5_files):
#       f = h5py.File('/home/madhu/datasets/download/test_del_analysis/new.fast5', 'r')
      f = h5py.File(a_fast5_file, 'r')
      if 'Analyses' in list(f.keys()):
         print('This file contains OLD Analysis: ', a_fast5_file)
         found_old_analysis = True
      else:
         print('This file is fine: ', a_fast5_file)

#    basefolder = '/home/madhu/datasets/download/guppy_v5.0.11/workspace/m6A_RNA_modification_native_RNA_seq/GSM3528749/data/extracted_data/fast5'
   
#    print('This is the argv: ', sys.argv)
#    print('\n\nThe length of argv: ', len(sys.argv))
   
#    print('This is the basefolder: ', basefolder)

#    if int(sys.argv[2])<1:
#       f5fs = glob.glob(os.path.join(basefolder, "*.fast5") )
#    else:
#       f5fs = glob.glob(os.path.join(basefolder, '/'.join(['*' for _ in range(int(sys.argv[2]))]) )) 
   
#    pmanager = multiprocessing.Manager();
#    handlers = []
#    h5files_Q = pmanager.Queue();

#    for h5f in f5fs:
#       h5files_Q.put( h5f )

#    share_var = (int(sys.argv[2]), h5files_Q )

#    print("{}".format( f5fs ))
#    #sys.exit(1)

#    for hid in range( int(sys.argv[3] )):
#       p = multiprocessing.Process(target=del_analysis, args=share_var);
#       p.start();
#       handlers.append(p);
   
#    while any(p.is_alive() for p in handlers):
#       try:
#          pass;
#       except:
#          time.sleep(1);
#          continue;
   if found_old_analysis == True: 
      print("\n\nXXXXXXX There are files with OLD Analysis. XXXXXXX\n\n")
   else:
      print('\n\nAll Good.\n\n')

