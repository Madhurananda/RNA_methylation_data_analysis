import os,sys
import glob
import tarfile
from tqdm import tqdm

# calculate file size in KB, MB, GB
def convert_bytes(size):
    """ Convert bytes to KB, or MB or GB"""
    for x in ['bytes', 'KB', 'MB', 'GB', 'TB']:
        if size < 1024.0:
            return "%3.1f %s" % (size, x)
        size /= 1024.0

def track_progress(members):
  for member in members:
     # this will be the current file being extracted
     yield member
    
        
if __name__=='__main__':
   if len(sys.argv)<3:
      print("Usage: {} {} {}".format(sys.argv[0], "basefolder", "output_folder"))
      print("Example:")
#       print("       {} {} {}".format(sys.argv[0], "Need to Update"))
#       print("       {} {} {}".format(sys.argv[0], "Curlcake_m6a", "1368282"))
#       print("       {} {} {}".format(sys.argv[0], "Curlcake_non", "1204670"))
#       print("       {} {} {}".format(sys.argv[0], "Curlcake_non", "1368288"))
      sys.exit(1)

   basefolder = sys.argv[1]
   output_folder = sys.argv[2]
   
   if basefolder.endswith('/'):
      print('')
   else:
      basefolder+='/'
    
   print('\n\nbasefolder:', basefolder)
   print('output_folder: ', output_folder)
   
#    data_dir = '../data/S_datasets_3/ncbi_dataset/data/'



   compressed_files = list(glob.iglob(basefolder+'*.tar.*', recursive=False))
#    compressed_files = list(glob.iglob(basefolder+'original_data/*.tar.*', recursive=False))
   
   if len(compressed_files) != len(set(compressed_files)):
      print('\n\n#### There are dupliacted tar files in the original data. Exiting ... #####\n\n')
      sys.exit(1)

   print('The number of compressed tar files: ', len(compressed_files))
#    print('The list of compressed files: ', compressed_files)

   for a_com_file in tqdm(compressed_files):
      print('\n\n')
#       a_com_file = compressed_files[0]
#       output_folder = basefolder+'extracted_data/'
    #    extract_dir = output_folder+a_com_file.split('.tar')[0]
      extract_dir = output_folder+a_com_file.split('.tar')[0].split('/')[-1]

      if not os.path.isdir(extract_dir):
         print('**\nThe Folder has been created: ', extract_dir, '**\n')
         os.system("mkdir {}".format( extract_dir ));

         print('Uncompressed file: ', a_com_file, ' of size: ', convert_bytes(os.path.getsize(a_com_file)) )
         print('destination: ', extract_dir)

         file = tarfile.open(a_com_file)
         file.extractall(extract_dir)
         print('\n\n')
      else:
         print('\n\n\n#### The directory already exists: ', extract_dir, '## Skipping ... ###\n\n')



