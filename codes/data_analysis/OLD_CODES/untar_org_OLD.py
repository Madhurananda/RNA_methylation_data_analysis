import os,sys
sys.path.insert(0, '/home/madhu/datasets/codes')

import glob
import tarfile
import gzip
from tqdm import tqdm

from datetime import datetime
from main import check_endWith, check_create_dir, check_conda_env

# calculate file size in KB, MB, GB
def convert_bytes(size):
    """ Convert bytes to KB, or MB or GB"""
    for x in ['bytes', 'KB', 'MB', 'GB', 'TB']:
        if size < 1024.0:
            return "%3.1f %s" % (size, x)
        size /= 1024.0

    
        
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
   
   basefolder = sys.argv[1]
   output_folder = sys.argv[2]
   zip_type = int(sys.argv[3])
   
   
   basefolder = check_endWith(basefolder)

   now = datetime.now()
   current_time = now.strftime("%Y/%m/%d at %H:%M:%S")
   print('\n\nThe script starting at: ' + str(current_time), ' \n\n' )
   
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
   
   for a_com_file in tqdm(compressed_files):
      print('\n\n')
      split_find_express = find_express.replace('*', '')
      extract_dir = output_folder+a_com_file.split(split_find_express)[0].split('/')[-1]
      
      ## If the folder already exists, just skip the file ... 
      if not os.path.isdir(extract_dir):
         print('**\nThe Folder has been created: ', extract_dir, '**\n')
         os.system("mkdir {}".format( extract_dir ));
         
         print('Uncompressed file: ', a_com_file, ' of size: ', convert_bytes(os.path.getsize(a_com_file)) )
         print('\nDestination: ', extract_dir)
         
         if zip_type == 1:
            out_file = extract_dir+"/"+a_com_file.split(split_find_express)[0].split('/')[-1]+'.fast5'
#             print('The output file is: ', out_file)
            input = gzip.GzipFile(a_com_file, 'rb')
            s = input.read()
            input.close()
#             output = open(extract_dir+"/"+a_com_file.replace(split_find_express, ''), 'wb')
            output = open(out_file, 'wb')
            output.write(s)
            output.close()
         elif zip_type == 2 or zip_type == 3 or zip_type == 4:
            file = tarfile.open(a_com_file)
            file.extractall(extract_dir)
         print('\n\n')
      else:
         print('\n\n\n#### The directory already exists: ', extract_dir, '## Skipping ... ###\n\n\n')
   
   now = datetime.now()
   current_time = now.strftime("%Y/%m/%d at %H:%M:%S")
   print('\n\nThe script completed at: ' + str(current_time), ' \n\n' )


