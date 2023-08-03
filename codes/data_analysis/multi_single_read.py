import os,sys
sys.path.insert(0, '/home/madhu/datasets/codes')
from main import check_endWith, check_create_dir, check_conda_env

from tqdm import tqdm

import numpy as np
from datetime import datetime
from multiprocessing.pool import ThreadPool, Pool
from multiprocessing import cpu_count

def do_multi_single(args):
   t0 = datetime.now()
   a_dir_part, basefolder, output_folder = args[0], args[1], args[2]
   try: 
      ex_dir_part = a_dir_part.split(basefolder)[-1]
      print('The ex_dir_part part is: ', ex_dir_part)
      save_dir = output_folder + ex_dir_part
      print('The save_dir part is: ', save_dir)

      check_create_dir(save_dir)
      multi_single_cmd = "multi_to_single_fast5 --input {} --save_path {} --recursive".format( a_dir_part, save_dir )
      print('The multiple to single command is: ', multi_single_cmd, ' for: ', a_dir_part)
      sys.stdout.flush()
      stream = os.popen(multi_single_cmd)
      output = stream.read()
      # os.system( multi_single_cmd )
      print('\n\nThe output: ------------------------------------------------ \n', (output), '\n\n')
      sys.stdout.flush()

   except Exception as e:
      print('a_dir_part: ', a_dir_part)
      print('Exception:', e)

   sys.stdout.flush()
   return(a_dir_part, datetime.now() - t0)


if __name__=='__main__':

   check_conda_env('ACONDA')

   if len(sys.argv)<3:
      print('Usage: {} {} {}'.format('multi_single_read.py', 'basefolder', 'output_folder'))
      print("Example:")
      print('       {} {} {} {}'.format('/home/madhu/datasets/scripts/chris_scripts/untar_org.py', "/home/madhu/datasets/download/m6A_RNA_modification_native_RNA_seq/GSM3528751/original_data/extracted_data/fast5/", "/home/madhu/datasets/download/m6A_RNA_modification_native_RNA_seq/GSM3528751/splitted_extracted_data/", '2'))
      print('\n\nThe basefolders should contain the zip files. There are three types of zipped files: 1: *.fast5.gz, 2: *.fast5.tar, 3: *.tar.gz* \n\n')
#       print("       {} {} {}".format(sys.argv[0], "Curlcake_m6a", "1368282"))
#       print("       {} {} {}".format(sys.argv[0], "Curlcake_non", "1204670"))
#       print("       {} {} {}".format(sys.argv[0], "Curlcake_non", "1368288"))
      sys.exit(1)

   startTime = datetime.now()
   current_time = startTime.strftime("%Y/%m/%d at %H:%M:%S")
   print('\n\nThe script starting at: ' + str(current_time), ' \n\n' )

   basefolder = sys.argv[1]
   output_folder = sys.argv[2]

   basefolder = check_endWith(basefolder)
   output_folder = check_endWith(output_folder)

   print('\n\nbasefolder:', basefolder)
   print('output_folder: ', output_folder)
   
   fast5files = [os.path.join(root, name)
             for root, dirs, files in os.walk(basefolder)
             for name in files
             if name.endswith((".fast5"))]
   
   print('The number of FAST5 files: ', len(fast5files), ' inside: ', basefolder)

   list_dir_part = []
   for a_fast5files in tqdm(fast5files):
#       print('The number of compressed zipped files: ', len(fast5files))
#       print('The entire file is : ', a_fast5files)
      a_fast5files_part = a_fast5files.split('/')[-1]
#       print('The ONLY file is: ', a_fast5files_part)
      dir_part = a_fast5files.split(a_fast5files_part)[0]
#       print('The dir part is: ', dir_part)
      list_dir_part.append(dir_part)
   
   unique_list_dir_part = list(np.unique(np.array(list_dir_part)))
   # print('The original len: ', len(list_dir_part))
   # print('The unique lan: ', len(unique_list_dir_part))
   
   inputs = zip(unique_list_dir_part, [basefolder]*len(unique_list_dir_part), [output_folder]*len(unique_list_dir_part))

   ## Use the maximum CPU possible 
   cpus = cpu_count()
   n_jobs = min(cpus, len(unique_list_dir_part))

   # n_jobs = 20

   print('n_jobs: ', n_jobs)
   sys.stdout.flush()
   ## I am using multiple Processors as it is faster than threads. 
   results = tqdm(Pool(n_jobs).imap_unordered(do_multi_single, inputs), total=len(unique_list_dir_part))
   # results = ThreadPool(n_jobs).imap_unordered(unzip, inputs)
   for result in results:
      # print('time (s):', result)
      print('File: ', result[0], 'time :', result[1])


   # print('\n\n')
   executionTime = (datetime.now() - startTime)
   current_time = datetime.now().strftime("%Y/%m/%d at %H:%M:%S")
   print('\n\nThe script completed at: ' + str(current_time))
   print('Execution time: ' + str(executionTime), ' \n\n')

