import os,sys
import glob
import tarfile
import gzip
from tqdm import tqdm

# calculate file size in KB, MB, GB
def convert_bytes(size):
    """ Convert bytes to KB, or MB or GB"""
    for x in ['bytes', 'KB', 'MB', 'GB', 'TB']:
        if size < 1024.0:
            return "%3.1f %s" % (size, x)
        size /= 1024.0

    
        
if __name__=='__main__':
   if len(sys.argv)<3:
      print('Usage: {} {} {} {}'.format('untar_org.py', 'basefolder', 'output_folder', 'zip_type'))
      print("Example:")
      print('       {} {} {} {}'.format('/home/madhu/datasets/scripts/chris_scripts/untar_org.py', "/home/madhu/datasets/download/m6A_RNA_modification_native_RNA_seq/GSM3528751/original_data/extracted_data/fast5/", "/home/madhu/datasets/download/m6A_RNA_modification_native_RNA_seq/GSM3528751/splitted_extracted_data/", '2'))
      print('\n\nThe basefolders should contain the zip files. There are three types of zipped files: 1: *.fast5.gz, 2: *.fast5.tar, 3: *.tar.gz* \n\n')
#       print("       {} {} {}".format(sys.argv[0], "Curlcake_m6a", "1368282"))
#       print("       {} {} {}".format(sys.argv[0], "Curlcake_non", "1204670"))
#       print("       {} {} {}".format(sys.argv[0], "Curlcake_non", "1368288"))
      sys.exit(1)
   
   basefolder = sys.argv[1]
   output_folder = sys.argv[2]
#    zip_type = int(sys.argv[3])
   
   if basefolder.endswith('/'):
      print('')
   else:
      basefolder+='/'
   
   print('\n\nbasefolder:', basefolder)
   print('output_folder: ', output_folder)
   
   fast5files = [os.path.join(root, name)
             for root, dirs, files in os.walk(basefolder)
             for name in files
             if name.endswith((".fast5"))]

   for a_fast5files in tqdm(fast5files[0:5]):
      print('The number of compressed zipped files: ', len(fast5files))
#    print('The files are: ', fast5files)
#       a_fast5files = fast5files[0]
      print('The entire file is : ', a_fast5files)
      a_fast5files_part = a_fast5files.split('/')[-1]
      print('The ONLY file is: ', a_fast5files_part)
      ex_dir_part = a_fast5files.split(basefolder)[-1]
      dir_part = a_fast5files.split(a_fast5files_part)[0]
      print('The dir part is: ', dir_part)
      print('The ex_dir_part part is: ', ex_dir_part)
      ext_file = output_folder + ex_dir_part
      print('The ext_file is: ', ext_file)
      save_dir = ext_file.split(a_fast5files_part)[0]
      print('The save_dir part is: ', save_dir)

      if not os.path.isdir(save_dir):
         os.makedirs(save_dir, exist_ok=True)
         print('**\nThe Folder has been created: ', save_dir, '**\n')
#       os.system("mkdir {}".format( save_dir ));
      else:
         print('The folder already exists: ', save_dir)
      
# #    os.system( "tar -xf {} -C {}".format( _f, basefolder+'/'+basetar+'/'+_fid ) )
#       os.system( "multi_to_single_fast5 --input {} --save_path {} --recursive".format( a_fast5files, save_dir ) )
      os.system( "multi_to_single_fast5 --input {} --save_path {}".format( a_fast5files, save_dir ) )




#    ## There are three types of zipped files: 1: *.fast5.gz, 2: *.fast5.tar, 3: *.tar.gz* 
#    if zip_type == 1:
#       find_express = '*.fast5.gz'
#    elif zip_type == 2:
#       find_express = '*.fast5.tar'
#    elif zip_type == 3:
#       find_express = '*.tar.gz*'
   
#    compressed_files = list(glob.iglob(basefolder+find_express, recursive=True))
#    compressed_files = list(glob.iglob(basefolder+'original_data/*.tar.*', recursive=False))
   
#    if len(compressed_files) != len(set(compressed_files)):
#       print('\n\n#### There are dupliacted tar files in the original data. Exiting ... #####\n\n')
#       sys.exit(1)
   
#    print('The number of compressed zipped files: ', len(compressed_files))
# #    print('The list of compressed files: ', compressed_files)
   
#    for a_com_file in tqdm(compressed_files):
#       print('\n\n')
#       split_find_express = find_express.replace('*', '')
#       extract_dir = output_folder+a_com_file.split(split_find_express)[0].split('/')[-1]
      
#       ## If the folder already exists, just skip the file ... 
#       if not os.path.isdir(extract_dir):
#          print('**\nThe Folder has been created: ', extract_dir, '**\n')
#          os.system("mkdir {}".format( extract_dir ));
         
#          print('Uncompressed file: ', a_com_file, ' of size: ', convert_bytes(os.path.getsize(a_com_file)) )
#          print('\nDestination: ', extract_dir)
         
#          if zip_type == 1:
#             out_file = extract_dir+"/"+a_com_file.split(split_find_express)[0].split('/')[-1]+'.fast5'
# #             print('The output file is: ', out_file)
#             input = gzip.GzipFile(a_com_file, 'rb')
#             s = input.read()
#             input.close()
# #             output = open(extract_dir+"/"+a_com_file.replace(split_find_express, ''), 'wb')
#             output = open(out_file, 'wb')
#             output.write(s)
#             output.close()
#          elif zip_type == 2 or zip_type == 3:
#             file = tarfile.open(a_com_file)
#             file.extractall(extract_dir)
#          print('\n\n')
#       else:
#          print('\n\n\n#### The directory already exists: ', extract_dir, '## Skipping ... ###\n\n\n')


   print('\n\n')
