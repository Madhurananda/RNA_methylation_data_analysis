import os
import sys
sys.path.insert(0, '/home/madhu/datasets/codes')

from main import check_endWith, check_create_dir, check_conda_env
from datetime import datetime
from glob import glob
from tqdm import tqdm
from multiprocessing.pool import ThreadPool, Pool
from multiprocessing import cpu_count


'''
This script generates the commands necessary for 'generate_tfrecord_bp.py'

input: python do_generate_tfrec_bp.py parent_basecall_dir parent_analysis_dir thread, single_mult, len_thr, map_qual_thr 

python /home/madhu/datasets/codes/ML_codes/madhu_codes/do_generate_tfrec_bp.py /home/madhu/datasets/basecall/guppy_v5.0.11/p3_m6A_RNA_modification_native_RNA_seq /home/madhu/datasets/analysis/p3_m6A_RNA_modification_native_RNA_seq 10 2 200 20 RNA > /home/madhu/LOGs/gen_TFRecbp.05m.p3_m6A_RNA_modification_native_RNA_seq.txt

'''

def do_TFRec(args):
    t0 = datetime.now()
    a_dir, basecalled_folder, analysis_folder  = args[0], args[1], args[2]
    try:
        current_basecall_folder = basecalled_folder + a_dir + '/'
        current_analysis_folder = analysis_folder + a_dir + '/'
        # current_basecall_folder = check_endWith(current_basecall_folder)
        # current_analysis_folder = check_endWith(current_analysis_folder)
        current_TFRec_dir = current_analysis_folder+'TFRec_BP/'
        BAM_file = current_analysis_folder+a_dir+'.bam'
        # print('current_basecall_folder: ', current_basecall_folder)
        # print('current_analysis_folder: ', current_analysis_folder)
        # print('current_TFRec_dir: ', current_TFRec_dir)
        # print('BAM_file: ', BAM_file)

        ## skip the folders with TFRec as a subfolder 
        if os.path.exists(current_TFRec_dir):
            print('\n\n****The TFRecs have already been generated for: ', current_analysis_folder, ' ****')
        else: 
            TFRec_cmd = 'python /home/madhu/datasets/codes/ML_codes/madhu_codes/generate_tfrecord_bp.py {} {} {} {} "workspace/*/" {} {} {} RNA'.format(current_basecall_folder, current_TFRec_dir, thread, single_mult, BAM_file, len_thr, map_qual_thr)

            print('\n\nTF-Rec BP command: ', TFRec_cmd, '\n\n')

            ## Execute the TFRec part here ... 
            os.system(TFRec_cmd)
            print('\n\n'+TFRec_cmd + ' has been successfully executed \n\n')

    except Exception as e:
        print('a_dir: ', a_dir)
        print('Exception:', e)
    
    sys.stdout.flush()
    return(a_dir, datetime.now() - t0)


if __name__=='__main__':
    if len(sys.argv)<2:
        print("Usage: {} {}".format(sys.argv[0], "basefolder"))
        print("Example:")
        print("       python {} {} {} {}".format(sys.argv[0], "/home/madhu/datasets/download/m6A_RNA_modification_native_RNA_seq/GSM3528751/extracted_data"))
#           print("       python {} {} {} {}".format(sys.argv[0], "Curlcake_non"))
        sys.exit(1)
    
    check_conda_env('ACONDA')

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

    analysis_folder = sys.argv[2]
    analysis_folder = check_endWith(analysis_folder)

    thread = sys.argv[3]
    single_mult = sys.argv[4]
    len_thr = sys.argv[5]
    map_qual_thr = sys.argv[6]

    ## Get all the subdirectoires of the base folder 
    dir_list = next(os.walk(basecalled_folder))[1]
    # print(intial_dir_list)
    

    inputs = zip(dir_list, [basecalled_folder]*len(dir_list), [analysis_folder]*len(dir_list))

    ## Use the maximum CPU possible 
    cpus = cpu_count()
    n_jobs = min(cpus, len(dir_list))

    # n_jobs = 20

    print('n_jobs: ', n_jobs)
    sys.stdout.flush()
    ## I am using multiple Processors as it is faster than threads. 
    results = tqdm(Pool(n_jobs).imap_unordered(do_TFRec, inputs), total=len(dir_list))
    # results = ThreadPool(n_jobs).imap_unordered(unzip, inputs)
    for result in results:
        # print('time (s):', result)
        print('File: ', result[0], 'time :', result[1])
    
    
    executionTime = (datetime.now() - startTime)
    current_time = datetime.now().strftime("%Y/%m/%d at %H:%M:%S")
    print('\n\nThe script completed at: ' + str(current_time))
    print('Execution time: ' + str(executionTime), ' \n\n')


