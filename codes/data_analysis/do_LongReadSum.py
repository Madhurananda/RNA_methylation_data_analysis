import os
import sys
sys.path.insert(0, '/home/madhu/work/codes')

from main import check_endWith, check_create_dir, check_env
from tqdm import tqdm
from datetime import datetime
from multiprocessing.pool import ThreadPool, Pool
from multiprocessing import cpu_count

'''
This script will generate LongReadSum analysis for each dataset. 
This script needs to be run from lrst_py39 environment. This is why it is a seperate script. 
This script checks the BAM inside the given folder and saves the output .... 
'''

def do_LongRead_analysis(args):
    t0 = datetime.now()
    BAM_file = args[0]
    try:
        output_dir = BAM_file.replace(BAM_file.split('/')[ len(BAM_file.split('/'))-1], 'long_output/') 
        # output_dir = base_folder + BAM_file.split('/')[ len(BAM_file.split('/'))-2]+'/long_output/'
        print('BAM_file: ', BAM_file)
        print('out_dir: ', output_dir)
        print('\n')
        # if check_create_dir(output_dir):
        longRD_cmd = 'python /home/madhu/work/codes/LongReadSum bam -i {} -o {}'.format(BAM_file, output_dir)
        # print('longRD_cmd: ', longRD_cmd)
        os.system(longRD_cmd)
        # else:
        #     print('Analysis has already been completed for: ', output_dir)
    except Exception as e:
        print('BAM_file: ', BAM_file)
        print('Exception:', e)
    
    sys.stdout.flush()
    pbar.update(1)
    return(BAM_file, datetime.now() - t0)


if __name__=='__main__':

    check_env('lrst_py39')

    startTime = datetime.now()
    current_time = startTime.strftime("%Y/%m/%d at %H:%M:%S")
    print('\n\nThe script starting at: ' + str(current_time), ' \n\n' )

    base_folder = sys.argv[1]
    base_folder = check_endWith(base_folder)

    # base_folder = '/home/madhu/datasets/analysis/p2_identification_m6A_m5C_RNA_modification'
    # base_folder = check_endWith(base_folder)

    all_BAM_files = [os.path.join(root, name)
                for root, dirs, files in os.walk(base_folder)
                for name in files
                if name.endswith((".bam"))]

    # print('all_BAM_files: ', all_BAM_files)

    inputs = zip(all_BAM_files)
    pbar = tqdm(total=len(all_BAM_files))

    ## Use the maximum CPU possible 
    cpus = cpu_count()
    n_jobs = min(cpus, len(all_BAM_files))

    # n_jobs = 20

    print('n_jobs: ', n_jobs)
    ## I am using multiple Processors as it is faster than threads. 
    results = tqdm(Pool(n_jobs).imap_unordered(do_LongRead_analysis, inputs), total=len(all_BAM_files))
    # results = ThreadPool(n_jobs).imap_unordered(unzip, inputs)
    for result in results:
        # print('time (s):', result)
        print('File: ', result[0], 'time :', result[1])
    


    # for BAM_file in tqdm(all_BAM_files):
    #     output_dir = BAM_file.replace(BAM_file.split('/')[ len(BAM_file.split('/'))-1], 'long_output/') 
    #     # output_dir = base_folder + BAM_file.split('/')[ len(BAM_file.split('/'))-2]+'/long_output/'
    #     print('BAM_file: ', BAM_file)
    #     print('out_dir: ', output_dir)
    #     print('\n')
    #     if check_create_dir(output_dir):
    #         longRD_cmd = 'python /home/madhu/datasets/codes/LongReadSum bam -i {} -o {}'.format(BAM_file, output_dir)
    #         # print('longRD_cmd: ', longRD_cmd)
    #         os.system(longRD_cmd)
    #     else:
    #         print('Analysis has already been completed for: ', output_dir)


    executionTime = (datetime.now() - startTime)
    current_time = datetime.now().strftime("%Y/%m/%d at %H:%M:%S")
    print('\n\nThe script completed at: ' + str(current_time))
    print('Execution time: ' + str(executionTime), ' \n\n')
