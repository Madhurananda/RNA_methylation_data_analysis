import os
import sys
import pickle

def check_endWith(basecalled_folder):
    if not basecalled_folder.endswith('/'):
        basecalled_folder+='/'
    return basecalled_folder

def check_create_dir(folder):
    if_create = False
    if not os.path.isdir(folder):
        os.makedirs(folder, exist_ok=True)
        if_create = True
        print('**\nThe Folder has been created: ', folder, '**\n')
    else:
        print('The folder already exists: ', folder)
    return if_create

def check_env(env):
    if not os.environ['CONDA_DEFAULT_ENV'] == env:
        print('\n\n*** Please use:    conda activate '+env+'    and run the script again.***\n\n')
        sys.exit(1)

# This function writes content to a file. 
def write_file(file_name, to_write):
    # replaces the contents 
    file = open(file_name,"w")#write mode 
    file.write(to_write) 
    file.close() 

# This function appends content to a file. 
def append_file(file_name, to_write):
    # Append-adds at last 
    file = open(file_name,"a")#append mode 
    file.write(to_write) 
    file.close() 

# calculate file size in KB, MB, GB
def convert_bytes(size):
    """ Convert bytes to KB, or MB or GB"""
    for x in ['bytes', 'KB', 'MB', 'GB', 'TB']:
        if size < 1024.0:
            return "%3.1f %s" % (size, x)
        size /= 1024.0

def count_lines(a_fastQ_file):
    stream = os.popen('cat {} | wc -l'.format(a_fastQ_file))
    output = stream.read()
#    print('line no: ', output, ' for file: ', a_fastQ_file)
    return output

def del_dir_content_files(dir):
    delete_dir_content_cmd = 'rm '+dir+'*'
    os.system(delete_dir_content_cmd)
    print('The files are deleted inside: ', dir, ' using the command: ', delete_dir_content_cmd)
    sys.stdout.flush()

def read_file(file):
    f = open(file, "r")
    return f.read()



def save_file(file_save_name, file_data):
    # with open('/home/madhu/work/comp_tools/paper_6/download_tools/EpiNano/df.pickle', 'wb') as handle:
    with open(file_save_name, 'wb') as handle:
        pickle.dump(file_data, handle, protocol=pickle.HIGHEST_PROTOCOL)

def load_file(file_load_name):
    # with open('/home/madhu/work/comp_tools/paper_6/test_paper_6/m6A_detection/pd.pkl', 'rb') as file:
    with open(file_load_name, 'rb') as file:
        file_data = pickle.load(file)
    
    return file_data

