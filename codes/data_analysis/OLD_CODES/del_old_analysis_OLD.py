
import os;
import sys;
import string;
import glob;
import time
import copy

from tqdm import tqdm

from datetime import datetime

import h5py
import numpy as np
import multiprocessing

from collections import defaultdict

import tempfile
import subprocess

import re;


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

def del_analysis( h5files_Q):
   while not h5files_Q.empty():
      try:
#          print('\n\n\nThe queue is not EMPTY.')
         f5t = h5files_Q.get(block=False)
#          print('\n\n\nThe length of the f5files, ', len(f5files))
#          print('\n\n\nOne f5files, ', f5files[0])
      except:
#          print('\n\n\nThe queue is EMPTY.')
         break;

      try:
         with h5py.File(f5t, 'r+') as mf5:
            if '/Analyses' in mf5:
               del mf5['/Analyses']
               mf5.flush();
               mf5.close();
               print('This file contained OLD analyses, which has been removed: ', f5t)
            else:
               print('\n\nThis file does not contain OLD analyses: ', f5t)
      except:
          mf5.close();
          print("Failed to delete analysis for {}".format(f5t))
    


if __name__=='__main__':
   if len(sys.argv)<3:
      print("Usage: {} {} {}".format(sys.argv[0], "basefolder", "threads"))
      print("Example:")
      print("       python {} {} {} {}".format(sys.argv[0], "Curlcake_m6a", 2, 7))
      print("       python {} {} {} {}".format(sys.argv[0], "Curlcake_non", 2, 7))
      sys.exit(1)
   
   basefolder = sys.argv[1]
    
   startTime = datetime.now()
   current_time = startTime.strftime("%Y/%m/%d at %H:%M:%S")
   print('\n\nThe script starting at: ' + str(current_time), ' \n\n' )
   
   
#    basefolder = '/home/madhu/datasets/download/guppy_v5.0.11/workspace/m6A_RNA_modification_native_RNA_seq/GSM3528749/data/extracted_data/fast5'
   
   
   f5fs = [os.path.join(root, name)
             for root, dirs, files in os.walk(basefolder)
             for name in files
             if name.endswith((".fast5"))]
   print('\n\nThe number of FST5 files are: ', len(f5fs), '\n\n')
#    print('A file-len is: ', len(f5fs[0]))

   pmanager = multiprocessing.Manager();
   handlers = []
   h5files_Q = pmanager.Queue();

   for h5f in f5fs:
      h5files_Q.put( h5f )



   for hid in range( int(sys.argv[2] )):
      p = multiprocessing.Process(target=del_analysis, args=(h5files_Q,));
      p.start();
      handlers.append(p);
   
   while any(p.is_alive() for p in handlers):
      try:
         pass;
      except:
         time.sleep(1);
         continue;

   print("\nAll Done!")
   executionTime = (datetime.now() - startTime)
   current_time = datetime.now().strftime("%Y/%m/%d at %H:%M:%S")
   print('\n\nThe script completed at: ' + str(current_time))
   print('Execution time: ' + str(executionTime), ' \n\n')
