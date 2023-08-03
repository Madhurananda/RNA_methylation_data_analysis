
import os;
import sys;
sys.path.insert(0, '/home/madhu/datasets/codes')

import string;
import glob;
import time
import copy

from tqdm import tqdm

import h5py
import numpy as np
from multiprocessing.pool import ThreadPool, Pool
from multiprocessing import cpu_count

from collections import defaultdict
from datetime import datetime
import tempfile
import subprocess

import re;
from main import check_endWith, check_create_dir, check_conda_env, write_file


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

def count_lines(a_fastQ_file):
   stream = os.popen('cat {} | wc -l'.format(a_fastQ_file))
   output = stream.read()
#    print('line no: ', output, ' for file: ', a_fastQ_file)
   return output

def do_check(args):
   t0 = datetime.now()

   try:
      to_wrtie_content = ''
      write_file_name = data_folder+'tombo_resquiggle_check_output.txt'

      ## Step 1: Read all the fast5 files in the basefolder 
      all_fast5_files = [os.path.join(root, name)
            for root, dirs, files in os.walk(data_folder)
            for name in files
            if name.endswith((".fast5"))]
      print('Len of all_fast5_files: ', len(all_fast5_files))
      to_wrtie_content += ('\nLen of all_fast5_files: '+str(len(all_fast5_files)) )
      sys.stdout.flush()

      ## Step 2: Calculate Tombo performance by analysing the FAST5 files 
      n_tombo = 0
      n_Ntombo =0
      for a_fast5_file in tqdm(all_fast5_files):
         f5r = h5py.File(a_fast5_file, 'r')
         if '/Analyses/RawGenomeCorrected_000/BaseCalled_template/Events' in f5r:
      #          print('This file is done with Tombo: ', a_fast5_file)
               n_tombo += 1
         else:
      #          print('This file is NOT done with Tombo: ', a_fast5_file)
               n_Ntombo += 1
      
      
      print('\n\nThe no of files with Tombo: ', n_tombo)
      print('\n\nThe no of files with NO Tombo: ', n_Ntombo)
      
      print('\n\n\n', a_dataset, ': The percentage of files for which tombo done: ', round(((n_tombo/(n_tombo+n_Ntombo))*100), 2), '%\n\n')
      
      to_wrtie_content += ('\nThe no of files with Tombo: '+str(n_tombo) )
      to_wrtie_content += ('\nThe no of files with NO Tombo: '+str(n_Ntombo) )
      to_wrtie_content += ('\n\n'+a_dataset+': The percentage of files for which tombo done:  '+str(round(((n_tombo/(n_tombo+n_Ntombo))*100), 2)) +'%\n\n')
      sys.stdout.flush()

      # Step 3:  Get all the FASTQ files
      all_fastQ_files = [os.path.join(root, name)
                  for root, dirs, files in os.walk(data_folder)
                  for name in files
                  if name.endswith((".fastq"))]
      print('Len of all_fastQ_files: ', len(all_fastQ_files))
      to_wrtie_content += ('\nLen of all_fastQ_files: '+str(len(all_fastQ_files)) )
      sys.stdout.flush()

      # Step 4:  Get the number of lines from each of those FASTQ files
      c_pass = 0
      c_fail = 0
      for a_fastQ_file in tqdm(all_fastQ_files): 
      #       print('The fastQ file is: ', a_fastQ_file)
         if data_folder+'pass/' in a_fastQ_file:
      #          print('The PASS fastQ file is: ', a_fastQ_file)
               n_lines = int(count_lines(a_fastQ_file))
               c_pass += n_lines
         elif data_folder+'fail/' in a_fastQ_file:
      #          print('The FAIL fastQ file is: ', a_fastQ_file)
               n_lines = int(count_lines(a_fastQ_file))
               c_fail += n_lines
         else:
               print('\n\n***** There is something VERY wrong::', a_fastQ_file, ':: INVESTIGATE...******\n\n')
      sys.stdout.flush()


      n_read_pass = c_pass/4
      n_read_fail = c_fail/4
      
      print('\nThe total pass line counts are: ', c_pass)
      print('\nThe total fail line counts are: ', c_fail)
      print('\n\n', a_dataset, ': pass/(pass+fail) READ percentage is: ', round((n_read_pass/(n_read_pass+n_read_fail))*100, 2), '%')

      to_wrtie_content += ('\nThe total pass line counts are: '+str(c_pass) )
      to_wrtie_content += ('\nThe total fail line counts are: '+str(c_fail) )
      to_wrtie_content += ('\n\n'+a_dataset+': pass/(pass+fail) READ percentage is: '+str(round((n_read_pass/(n_read_pass+n_read_fail))*100, 2)) + '%\n\n' )

      write_file(write_file_name, to_wrtie_content)

   except Exception as e:
      print('a_fast5_file: ', a_fast5_file)
      print('Exception:', e)
   
   sys.stdout.flush()
   return(a_dataset, datetime.now() - t0)



if __name__=='__main__':
   if len(sys.argv)<2:
      print("Usage: {} {}".format(sys.argv[0], "basefolder"))
      print("Example:")
      print("       python {} {} {} {}".format(sys.argv[0], "/home/madhu/datasets/download/m6A_RNA_modification_native_RNA_seq/GSM3528751/extracted_data"))
#       print("       python {} {} {} {}".format(sys.argv[0], "Curlcake_non"))
      sys.exit(1)
   
   basefolder = sys.argv[1]
   
   basefolder = check_endWith(basefolder)

   ## Step 1: Read all the fast5 files in the basefolder 
   
   all_fast5_files = [os.path.join(root, name)
             for root, dirs, files in os.walk(basefolder)
             for name in files
             if name.endswith((".fast5"))]
   print('Len of all_fast5_files: ', len(all_fast5_files))
   
   ## Step 2: Calculate Tombo performance by analysing the FAST5 files 
   n_tombo = 0
   n_Ntombo =0
   for a_fast5_file in tqdm(all_fast5_files):
      f5r = h5py.File(a_fast5_file, 'r')
      if '/Analyses/RawGenomeCorrected_000/BaseCalled_template/Events' in f5r:
#          print('This file is done with Tombo: ', a_fast5_file)
         n_tombo += 1
      else:
#          print('This file is NOT done with Tombo: ', a_fast5_file)
         n_Ntombo += 1
   
   
   print('\n\nThe no of files with Tombo: ', n_tombo)
   print('\n\nThe no of files with NO Tombo: ', n_Ntombo)
   
   print('\n\n\nThe percentage of tombo done: ', round(((n_tombo/(n_tombo+n_Ntombo))*100), 4), '%\n\n')
   
   # Step 3:  Get all the FASTQ files
   all_fastQ_files = [os.path.join(root, name)
             for root, dirs, files in os.walk(basefolder)
             for name in files
             if name.endswith((".fastq"))]
   print('Len of all_fastQ_files: ', len(all_fastQ_files))

   # Step 4:  Get the number of lines from each of those FASTQ files
   c_pass = 0
   c_fail = 0
   for a_fastQ_file in tqdm(all_fastQ_files): 
#       print('The fastQ file is: ', a_fastQ_file)
      if basefolder+'pass/' in a_fastQ_file:
#          print('The PASS fastQ file is: ', a_fastQ_file)
         n_lines = int(count_lines(a_fastQ_file))
         c_pass += n_lines
      elif basefolder+'fail/' in a_fastQ_file:
#          print('The FAIL fastQ file is: ', a_fastQ_file)
         n_lines = int(count_lines(a_fastQ_file))
         c_fail += n_lines
      else:
         print('\n\n***** There is something VERY wrong:::: INVESTIGATE...******\n\n')


   
   print('\nThe total pass line counts are: ', c_pass)
   print('\nThe total fail line counts are: ', c_fail)

   n_read_pass = c_pass/4
   n_read_fail = c_fail/4

   print('\n\npass/(pass+fail) READ percentage is: ', round((n_read_pass/(n_read_pass+n_read_fail))*100, 2), '%')

    







