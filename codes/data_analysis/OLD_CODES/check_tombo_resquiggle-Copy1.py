
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
Tombo added the following to the fast5 files .... 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

/Analyses/RawGenomeCorrected_000 Group
/Analyses/RawGenomeCorrected_000/BaseCalled_template Group
/Analyses/RawGenomeCorrected_000/BaseCalled_template/Alignment Group
/Analyses/RawGenomeCorrected_000/BaseCalled_template/Events Dataset {156}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

------ This one worked with Tombo ------

(base) madhu@qgenlabgpu:~$ h5ls -r /home/madhu/datasets/basecall/guppy_v5.0.11/m6A_RNA_modification_native_RNA_seq/GSM3528752/workspace/1368282-29/MinION_2_20181108_FAK35429_MN29580_sequencing_run_RNA081120182_38847_read_126919_ch_140_strand.fast5
/                        Group
/Analyses                Group
/Analyses/Basecall_1D_000 Group
/Analyses/Basecall_1D_000/BaseCalled_template Group
/Analyses/Basecall_1D_000/BaseCalled_template/Fastq Dataset {SCALAR}
/Analyses/Basecall_1D_000/BaseCalled_template/Move Dataset {5950}
/Analyses/Basecall_1D_000/BaseCalled_template/Trace Dataset {5950, 8}
/Analyses/Basecall_1D_000/Summary Group
/Analyses/Basecall_1D_000/Summary/basecall_1d_template Group
/Analyses/RawGenomeCorrected_000 Group
/Analyses/RawGenomeCorrected_000/BaseCalled_template Group
/Analyses/RawGenomeCorrected_000/BaseCalled_template/Alignment Group
/Analyses/RawGenomeCorrected_000/BaseCalled_template/Events Dataset {344}
/Analyses/Segmentation_000 Group
/Analyses/Segmentation_000/Summary Group
/Analyses/Segmentation_000/Summary/segmentation Group
/PreviousReadInfo        Group
/Raw                     Group
/Raw/Reads               Group
/Raw/Reads/Read_126919   Group
/Raw/Reads/Read_126919/Signal Dataset {59900/Inf}
/UniqueGlobalKey         Group
/UniqueGlobalKey/channel_id Group
/UniqueGlobalKey/context_tags Group
/UniqueGlobalKey/tracking_id Group



------ This one didn't worked with Tombo ------

(base) madhu@qgenlabgpu:~$ h5ls -r /home/madhu/datasets/basecall/guppy_v5.0.11/m6A_RNA_modification_native_RNA_seq/GSM3528752/workspace/1368282-29/MinION_2_20181108_FAK35429_MN29580_sequencing_run_RNA081120182_38847_read_40257_ch_18_strand.fast5
/                        Group
/Analyses                Group
/Analyses/Basecall_1D_000 Group
/Analyses/Basecall_1D_000/BaseCalled_template Group
/Analyses/Basecall_1D_000/BaseCalled_template/Fastq Dataset {SCALAR}
/Analyses/Basecall_1D_000/BaseCalled_template/Move Dataset {1934}
/Analyses/Basecall_1D_000/BaseCalled_template/Trace Dataset {1934, 8}
/Analyses/Basecall_1D_000/Summary Group
/Analyses/Basecall_1D_000/Summary/basecall_1d_template Group
/Analyses/RawGenomeCorrected_000 Group
/Analyses/RawGenomeCorrected_000/BaseCalled_template Group
/Analyses/Segmentation_000 Group
/Analyses/Segmentation_000/Summary Group
/Analyses/Segmentation_000/Summary/segmentation Group
/PreviousReadInfo        Group
/Raw                     Group
/Raw/Reads               Group
/Raw/Reads/Read_40257    Group
/Raw/Reads/Read_40257/Signal Dataset {25104/Inf}
/UniqueGlobalKey         Group
/UniqueGlobalKey/channel_id Group
/UniqueGlobalKey/context_tags Group
/UniqueGlobalKey/tracking_id Group


I am checking if all the fast5 files contain: "/Analyses/RawGenomeCorrected_000/BaseCalled_template/Events"

'''



if __name__=='__main__':
   if len(sys.argv)<2:
      print("Usage: {} {}".format(sys.argv[0], "basefolder"))
      print("Example:")
      print("       python {} {} {} {}".format(sys.argv[0], "/home/madhu/datasets/download/m6A_RNA_modification_native_RNA_seq/GSM3528751/extracted_data"))
#       print("       python {} {} {} {}".format(sys.argv[0], "Curlcake_non"))
      sys.exit(1)
   
   basefolder = sys.argv[1]
   
   if basefolder.endswith('/'):
      print('')
   else:
      basefolder+='/'
   

   ## Read all the fast5 files in the basefolder 
   
   all_fast5_files = [os.path.join(root, name)
             for root, dirs, files in os.walk(basefolder)
             for name in files
             if name.endswith((".fast5"))]
    
   
#    print('all_fast5_files: ', len(all_fast5_files))
   print('Len of all_fast5_files: ', len(all_fast5_files))
   
   n_tombo = 0
   n_Ntombo =0
   for a_fast5_file in tqdm(all_fast5_files):
      f5r = h5py.File(a_fast5_file, 'r')
      if '/Analyses/RawGenomeCorrected_000/BaseCalled_template/Events' in f5r:
         print('This file is done with Tombo: ', a_fast5_file)
         n_tombo += 1
      else:
         print('This file is NOT done with Tombo: ', a_fast5_file)
         n_Ntombo += 1

#       f = h5py.File(a_fast5_file, 'r')py
# #       print('f lists: ', list(f.keys()))
# #       print('f[Analyses] ', f['Analyses'].keys())
#       if 'Analyses' in list(f.keys()):
#          if 'RawGenomeCorrected_000' in list(f['Analyses'].keys()):
# #             print('This file is done with Tombo: ', a_fast5_file)
#             n_tombo += 1
#          else:
# #             print('This file is NOT done with Tombo: ', a_fast5_file)
#             n_Ntombo += 1

   print('\n\nThe no of files with Tombo: ', n_tombo)
   print('\n\nThe no of files with NO Tombo: ', n_Ntombo)

   print('\n\n\nThe percentage of tombo done: ', round(((n_tombo/(n_tombo+n_Ntombo))*100), 4), '%\n\n')



#    if found_old_analysis == True: 
#       print("\n\nXXXXXXX There are files with OLD Analysis. XXXXXXX\n\n")
#    else:
#       print('\n\nAll Good.\n\n')

