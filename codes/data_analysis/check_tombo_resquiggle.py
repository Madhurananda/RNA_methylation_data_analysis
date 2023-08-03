
import os;
import sys;
sys.path.insert(0, '/home/madhu/work/codes')

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
from main import check_endWith, check_create_dir, check_env, write_file, count_lines

'''
I have made this script to use multiprocessing in doing Tombo checking. I have used both Pool and ThreadPool here. 
I DO NOT need to run this script as it is very similar to the do_Tombo_resquiggle.py
Only time I might need to run it to investigate more about a dataset or if the Tombo ouptput is not enough. 
'''


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

def do_step_2(args):
   #  t0 = datetime.now()
    a_fast5_file = args[0]
    done = False
    try:
   #  for a_fast5_file in tqdm(all_fast5_files):
      f5r = h5py.File(a_fast5_file, 'r')
      if '/Analyses/RawGenomeCorrected_000/BaseCalled_template/Events' in f5r:
         done = True   
    except Exception as e:
      print('a_fast5_file: ', a_fast5_file)
      print('Exception:', e)
   
    sys.stdout.flush()
    return(done, a_fast5_file)


def do_step_4(args):
    a_fastQ_file, data_folder = args[0], args[1]
    in_pass = False
    try:
        n_lines = int(count_lines(a_fastQ_file))
        print('n_lines: ', n_lines)
        if data_folder+'pass/' in a_fastQ_file:
            in_pass = True
        elif data_folder+'fail/' in a_fastQ_file:
            in_pass = False
        else:
            print('\n\n***** There is something VERY wrong::', a_fastQ_file, ':: INVESTIGATE...******\n\n')
    except Exception as e:
        print('a_fastQ_file: ', a_fastQ_file)
        print('Exception:', e)
   
    sys.stdout.flush()
    print('in_pass: ', in_pass)
    return(in_pass, n_lines, a_fastQ_file)



def do_TomboCheck(args):
    t0 = datetime.now()
    a_dataset, basecalled_folder, n_CPU = args[0], args[1], args[2]
    # print('basecalled_folder: ', basecalled_folder)
    # print('n_CPU: ', n_CPU)
    try:
        # Make the same folder in the analysis directories  
        data_folder=basecalled_folder+a_dataset
        data_folder = check_endWith(data_folder)

      #   Tombo_cmd = 'tombo resquiggle {} {} --processes {} --overwrite --num-most-common-errors 5'.format(data_folder, ref_genome, n_CPU)
      #   print('\nTombo_cmd :', Tombo_cmd)
      #   os.system(Tombo_cmd)
      #   print('\n\n'+Tombo_cmd + ' has been successfully executed for: ', data_folder)
      #   sys.stdout.flush()

        '''
        Now, here also do the Tombo statistics ..., same as check_tombo_resquiggle.py .... 
        '''
        to_wrtie_content = ''
        write_file_name = data_folder+'tombo_resquiggle_check_output.txt'








        ## Step 1: Read all the fast5 files in the basefolder 
        all_fast5_files = [os.path.join(root, name)
             for root, dirs, files in os.walk(data_folder)
             for name in files
             if name.endswith((".fast5"))]
        print('Len of all_fast5_files: ', len(all_fast5_files), ' for: ', a_dataset)
        to_wrtie_content += '\nLen of all_fast5_files: '+str(len(all_fast5_files))+' for: '+str(a_dataset) 
        sys.stdout.flush()

        ## Step 2: Calculate Tombo performance by analysing the FAST5 files 
        n_tombo = 0
        n_Ntombo = 0
        
        inputs = zip(all_fast5_files)
        results = tqdm(Pool(int(n_CPU)).imap_unordered(do_step_2, inputs), total=len(all_fast5_files))
      #   print('\nIT  CAME HERE\n\n')      
      #   sys.stdout.flush()
        # results = ThreadPool(n_jobs).imap_unordered(unzip, inputs)
        for result in results:
            if result[0] == True:
               n_tombo += 1
            else:
               n_Ntombo += 1
        
        print('\n\n', a_dataset, ': The no of files with Tombo: ', n_tombo, ' for: ', a_dataset)
        print('\n\n', a_dataset, ': The no of files with NO Tombo: ', n_Ntombo, ' for: ', a_dataset)
        
        print('\n\n\n', a_dataset, ': The percentage of files for which tombo done: ', round(((n_tombo/(n_tombo+n_Ntombo))*100), 2), '%\n\n')
        
        to_wrtie_content += ('\n'+a_dataset+': The no of files with Tombo: '+str(n_tombo) )
        to_wrtie_content += ('\n'+a_dataset+': The no of files with NO Tombo: '+str(n_Ntombo) )
        to_wrtie_content += ('\n\n'+a_dataset+': The percentage of files for which tombo done:  '+str(round(((n_tombo/(n_tombo+n_Ntombo))*100), 2)) +'%\n\n')
        sys.stdout.flush()


        
        '''
        I was trying to calculate the pass/fail reads using multiprocessing, but it didn't work, so I am going back to the normal processing. 
        '''


      #   # Step 3:  Get all the FASTQ files
      #   all_fastQ_files = [os.path.join(root, name)
      #               for root, dirs, files in os.walk(data_folder)
      #               for name in files
      #               if name.endswith((".fastq"))]
      #   print('Len of all_fastQ_files: ', len(all_fastQ_files))
      #   to_wrtie_content += ('\nLen of all_fastQ_files: '+str(len(all_fastQ_files)) )
      #   sys.stdout.flush()

      #   # Step 4:  Get the number of lines from each of those FASTQ files
      #   c_pass = 0
      #   c_fail = 0

      #   inputs = zip(all_fastQ_files, [data_folder]*len(all_fastQ_files))
      #   results = tqdm(Pool(int(n_CPU)).imap_unordered(do_step_4, inputs), total=len(all_fastQ_files))
      # #   print('\n\nXXXXXXXX')
      # #   print(results)
      # #   print('XXXXXXXX\n\n')
      # #   print('results: ', len(results), '\n\n')
      #   if result in results:
      #       print('\n\nIT CAME HERE\n\n')
      #       print('result[0]: ', result[0])
      #       if result[0] == True:
      #           print('pass, n_line: ', result[1])
      #           c_pass += result[1]
      #       else:
      #           c_fail += result[1]
      #           print('fail, n_line: ', result[1])
        
      #   n_read_pass = c_pass/4
      #   n_read_fail = c_fail/4
        
      #   print('\n', a_dataset, ': The total pass line counts are: ', c_pass)
      #   print('\n', a_dataset, ': The total fail line counts are: ', c_fail)
      #   print('\n\n', a_dataset, ': pass/(pass+fail) READ percentage is: ', round((n_read_pass/(n_read_pass+n_read_fail))*100, 2), '%')

      #   to_wrtie_content += ('\n'+a_dataset+': The total pass line counts are: '+str(c_pass) )
      #   to_wrtie_content += ('\n'+a_dataset+': The total fail line counts are: '+str(c_fail) )
      #   to_wrtie_content += ('\n\n'+a_dataset+': pass/(pass+fail) READ percentage is: '+str(round((n_read_pass/(n_read_pass+n_read_fail))*100, 2)) + '%\n\n' )
        
      # #   write_file(write_file_name, to_wrtie_content)
      #   sys.stdout.flush()
        

        

        # Step 3:  Get all the FASTQ files
        all_fastQ_files = [os.path.join(root, name)
                    for root, dirs, files in os.walk(data_folder)
                    for name in files
                    if name.endswith((".fastq"))]
        print('Len of all_fastQ_files: ', len(all_fastQ_files), ' for: ', a_dataset)
        to_wrtie_content += '\nLen of all_fastQ_files: '+str(len(all_fastQ_files)) + ' for: ' +str(a_dataset)
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
        
        print('\n', a_dataset, ': The total pass reads are: ', str(int(n_read_pass)))
        print('\n', a_dataset, ': The total fail reads are: ', str(int(n_read_fail)))
        print('\n\n', a_dataset, ': pass/(pass+fail) READ percentage is: ', round((n_read_pass/(n_read_pass+n_read_fail))*100, 2), '%')

        to_wrtie_content += ('\n'+a_dataset+': The total pass line counts are: '+str(int(n_read_pass)) )
        to_wrtie_content += ('\n'+a_dataset+': The total fail line counts are: '+str(int(n_read_fail)) )
        to_wrtie_content += ('\n\n'+a_dataset+': pass/(pass+fail) READ percentage is: '+str(round((n_read_pass/(n_read_pass+n_read_fail))*100, 2)) + '%\n\n' )
        write_file(write_file_name, to_wrtie_content)
        sys.stdout.flush()

        # print('\n', a_dataset, ': The total pass line counts are: ', c_pass)
        # print('\n', a_dataset, ': The total fail line counts are: ', c_fail)
        # print('\n\n', a_dataset, ': pass/(pass+fail) READ percentage is: ', round((c_pass/(c_pass+c_fail))*100, 2), '%')

        # to_wrtie_content += ('\n'+a_dataset+': The total pass line counts are: '+str(c_pass) )
        # to_wrtie_content += ('\n'+a_dataset+': The total fail line counts are: '+str(c_fail) )
        # to_wrtie_content += ('\n\n'+a_dataset+': pass/(pass+fail) READ percentage is: '+str(round((c_pass/(c_pass+c_fail))*100, 2)) + '%\n\n' )
        # write_file(write_file_name, to_wrtie_content)
        # sys.stdout.flush()

    except Exception as e:
        print('a_dataset: ', a_dataset)
        print('Exception:', e)
    
    sys.stdout.flush()
    return(a_dataset, datetime.now() - t0)





if __name__=='__main__':
   if len(sys.argv)<2:
      print("Usage: {} {}".format(sys.argv[0], "basefolder"))
      print("Example:")
      print("       python {} {} {} {}".format(sys.argv[0], "/home/madhu/datasets/download/m6A_RNA_modification_native_RNA_seq/GSM3528751/extracted_data"))
#           print("       python {} {} {} {}".format(sys.argv[0], "Curlcake_non"))
      sys.exit(1)
   
   check_env('ont4')

   startTime = datetime.now()
   current_time = startTime.strftime("%Y/%m/%d at %H:%M:%S")
   print('\n\nThe script starting at: ' + str(current_time), ' \n\n' )

   ## Basecall folder should be the folder which has been used as the output of guppy. 
   # basecalled_folder = '/home/madhu/datasets/basecall/guppy_v5.0.11/p2_identification_m6A_m5C_RNA_modification/HeLa_all_fast5_M11A_KO'
   basecalled_folder = sys.argv[1]
   basecalled_folder = check_endWith(basecalled_folder)

   # dataset_list = sys.argv[2].split(',')
   # dataset_list = list(filter(None, dataset_list)) # remove empty strings 
   # # print(dataset_list)

   given_dataset_list = sys.argv[2]
   if given_dataset_list == 'ALL_IN_DIR':
      print('All the folders inside the given folder will be considered.')
      dataset_list = next(os.walk(basecalled_folder))[1]
   ## I can also make another option --exclude 
   else:
      dataset_list = given_dataset_list.split(',')
      dataset_list = list(filter(None, dataset_list)) # remove empty strings 
      # print(dataset_list)

   n_CPU = sys.argv[3]

   inputs = zip(dataset_list, [basecalled_folder]*len(dataset_list), [n_CPU]*len(dataset_list))

   ## Use the maximum CPU possible 
   cpus = cpu_count()
   n_jobs = min(cpus, len(dataset_list))

   # n_jobs = 20

   print('n_jobs: ', n_jobs)
   ## I am using multiple Processors as it is faster than threads. 
   # results = tqdm(Pool(n_jobs).imap_unordered(do_Tombo, inputs), total=len(dataset_list))
   results = tqdm(ThreadPool(n_jobs).imap_unordered(do_TomboCheck, inputs), total=len(dataset_list))
   for result in results:
      # print('time (s):', result)
      print('Dataset: ', result[0], 'time took :', result[1])
   
   executionTime = (datetime.now() - startTime)
   current_time = datetime.now().strftime("%Y/%m/%d at %H:%M:%S")
   print('\n\nThe script completed at: ' + str(current_time))
   print('Execution time: ' + str(executionTime), ' \n\n')


