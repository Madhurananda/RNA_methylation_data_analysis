import os,sys
sys.path.insert(0, '/home/madhu/work/codes')
from main import check_endWith, check_create_dir, check_env

# check_env('ACONDA')


import glob
import tarfile
import gzip
from tqdm import tqdm
# import time
from datetime import datetime
from multiprocessing.pool import ThreadPool, Pool
from multiprocessing import cpu_count


# calculate file size in KB, MB, GB
def convert_bytes(size):
    """ Convert bytes to KB, or MB or GB"""
    for x in ['bytes', 'KB', 'MB', 'GB', 'TB']:
        if size < 1024.0:
            return "%3.1f %s" % (size, x)
        size /= 1024.0

def unzip(args):
   t0 = datetime.now()
   a_com_file, output_folder = args[0], args[1]
   print('\n\n')
   
   # split_find_express = find_express.replace('*', '')
   extract_dir = output_folder+a_com_file.split('.')[0].split('/')[-1]
   
   # print('a_com_file: ', a_com_file)
   # print('extract_dir: ', extract_dir)

   try:
      ## If the folder already exists, just skip the file ... 
      if not os.path.isdir(extract_dir):
         print('**\nThe Folder has been created: ', extract_dir, '**\n')
         os.system("mkdir {}".format( extract_dir ));
         
         print('Uncompressed file: ', a_com_file, ' of size: ', convert_bytes(os.path.getsize(a_com_file)) )
         print('\nDestination: ', extract_dir, ' for: ', a_com_file)
         if '.tar' in a_com_file:
            # print('the zip types are 2 to 4')
            file = tarfile.open(a_com_file)
            file.extractall(extract_dir)
         elif ('.gz' in a_com_file):
            # print('this is type 1')
            # out_file = extract_dir+"/"+a_com_file.split(split_find_express)[0].split('/')[-1]+'.fast5'
            out_file = extract_dir+"/"+a_com_file.split('.gz')[0].split('/')[-1]
   #             print('The output file is: ', out_file)
            input = gzip.GzipFile(a_com_file, 'rb')
            s = input.read()
            input.close()
   #             output = open(extract_dir+"/"+a_com_file.replace(split_find_express, ''), 'wb')
            output = open(out_file, 'wb')
            output.write(s)
            output.close()
         print(a_com_file, ' has been unzipped at: ', str(datetime.now()))
         print('\n\n')
      else:
         print('\n\n\n#### The directory already exists: ', extract_dir, '## Skipping ... ###\n\n\n')

   except Exception as e:
        print('a_com_file: ', a_com_file)
        print('Exception:', e)
   
   sys.stdout.flush()
   return(a_com_file, datetime.now() - t0)


if __name__=='__main__':
   if len(sys.argv)<4:
      print('Usage: {} {} {} {}'.format('untar_org.py', 'basefolder', 'output_folder', 'zip_type'))
      print("Example:")
      print('       {} {} {} {}'.format('/home/madhu/datasets/codes/data_analysis/untar_org.py', "/home/madhu/datasets/download/m6A_RNA_modification_native_RNA_seq/GSM3528751/original_data/extracted_data/fast5/", "/home/madhu/datasets/download/m6A_RNA_modification_native_RNA_seq/GSM3528751/splitted_extracted_data/", '2'))
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
   zip_type = int(sys.argv[3])
   
   
   basefolder = check_endWith(basefolder)
   output_folder = check_endWith(output_folder)
   
   
   print('\n\nbasefolder:', basefolder)
   print('output_folder: ', output_folder)


   ## There are three types of zipped files: 1: *.fast5.gz, 2: *.fast5.tar, 3: *.tar.gz* 4: *.tar*
   if zip_type == 1:
      find_express = '*.fast5.gz'
   elif zip_type == 2:
      find_express = '*.fast5.tar'
   elif zip_type == 3:
      find_express = '*.tar.gz*'
   elif zip_type == 4:
      find_express = '*.tar*'

   compressed_files = list(glob.iglob(basefolder+find_express, recursive=False))
#    compressed_files = list(glob.iglob(basefolder+'original_data/*.tar.*', recursive=False))
   
   if len(compressed_files) != len(set(compressed_files)):
      print('\n\n#### There are dupliacted tar files in the original data. Exiting ... #####\n\n')
      sys.exit(1)
   
   print('The number of compressed zipped files: ', len(compressed_files))
#    print('The list of compressed files: ', compressed_files)
   
   # ## This is for testing only ..........
   # compressed_files = compressed_files[0:5]

   inputs = zip(compressed_files, [output_folder]*len(compressed_files))

   ## Use the maximum CPU possible 
   cpus = cpu_count()
   n_jobs = min(cpus, len(compressed_files))

   # n_jobs = 20

   print('n_jobs: ', n_jobs)
   sys.stdout.flush()
   ## I am using multiple Processors as it is faster than threads. 
   results = tqdm(Pool(n_jobs).imap_unordered(unzip, inputs), total=len(compressed_files))
   # results = ThreadPool(n_jobs).imap_unordered(unzip, inputs)
   for result in results:
      # print('time (s):', result)
      print('File: ', result[0], 'time :', result[1])


   executionTime = (datetime.now() - startTime)
   current_time = datetime.now().strftime("%Y/%m/%d at %H:%M:%S")
   print('\n\nThe script completed at: ' + str(current_time))
   print('Execution time: ' + str(executionTime), ' \n\n')


