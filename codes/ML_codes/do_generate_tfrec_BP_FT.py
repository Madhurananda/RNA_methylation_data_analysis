import os
import sys
sys.path.insert(0, '/home/madhu/work/codes')

from main import check_endWith, check_create_dir, check_conda_env
from datetime import datetime
from glob import glob
from tqdm import tqdm
from multiprocessing.pool import ThreadPool, Pool
from multiprocessing import cpu_count


'''
This script generates the commands necessary for 'generate_tfrecord_bp.py'

input: python do_generate_tfrec_bp.py parent_basecall_dir parent_analysis_dir thread, single_mult, len_thr, map_qual_thr 

python /home/madhu/datasets/codes/ML_codes/madhu_codes/do_generate_tfrec_bp.py /home/madhu/datasets/basecall/guppy_v5.0.11/p3_m6A_RNA_modification_native_RNA_seq /home/madhu/datasets/analysis/p3_m6A_RNA_modification_native_RNA_seq 10 2 200 20 RNA > /home/madhu/Logs/gen_TFRecbp.05m.p3_m6A_RNA_modification_native_RNA_seq.txt

'''



def do_TFRec(args):
    t0 = datetime.now()
    a_dir, basecalled_folder, analysis_folder  = args[0], args[1], args[2]
    try:
        current_basecall_folder = basecalled_folder + a_dir + '/'
        current_analysis_folder = analysis_folder + a_dir + '/'
        # current_basecall_folder = check_endWith(current_basecall_folder)
        # current_analysis_folder = check_endWith(current_analysis_folder)
        BAM_file = current_analysis_folder+a_dir+'.bam'

        if BP_FT_tougle == 'BP':
            current_TFRec_dir = current_analysis_folder+'TFRec_BP/'
            TFRec_cmd = 'python /home/madhu/work/codes/ML_codes/madhu_codes/generate_tfrecord_bp.py {} {} {} {} "workspace/*/" {} {} {} RNA'.format(current_basecall_folder, current_TFRec_dir, thread, single_mult, BAM_file, len_thr, map_qual_thr)
        elif BP_FT_tougle == 'FT':
            current_TFRec_dir = current_analysis_folder+'TFRec_FT/'
            TFRec_cmd = 'python /home/madhu/work/codes/ML_codes/madhu_codes/generate_tfrecord_ft.py {} {} {} {} "workspace/*/" {} {} {} RNA'.format(current_basecall_folder, current_TFRec_dir, thread, single_mult, BAM_file, len_thr, map_qual_thr)
        else:
            print('Please provide input: BP or FT')

        
        # print('current_basecall_folder: ', current_basecall_folder)
        # print('current_analysis_folder: ', current_analysis_folder)
        # print('current_TFRec_dir: ', current_TFRec_dir)
        # print('BAM_file: ', BAM_file)

        

        # ## skip the folders with TFRec as a subfolder 
        # if os.path.exists(current_TFRec_dir):
        #     print('\n\n****The TFRecs have already been generated for: ', current_analysis_folder, ' ****')
        # else: 
            
        print('\n\nTF-Rec command which will run now: ', TFRec_cmd, '\n\n')
        sys.stdout.flush()
        
        ## Execute the TFRec part here ... 
        stream = os.popen(TFRec_cmd)
        output = stream.read()
        print('\n\nThe output for: ', TFRec_cmd, ' is: \n', (output), '\n\n')

    except Exception as e:
        print('a_dir: ', a_dir)
        print('Exception:', e)
    
    sys.stdout.flush()
    return(a_dir, datetime.now() - t0)


if __name__=='__main__':
    if len(sys.argv)<2:
        print("Usage: {} {}".format(sys.argv[0], "basefolder"))
        print("Example:")
        print("       python {} {} {} {}".format(sys.argv[0], "/home/madhu/work/download/m6A_RNA_modification_native_RNA_seq/GSM3528751/extracted_data"))
#           print("       python {} {} {} {}".format(sys.argv[0], "Curlcake_non"))
        sys.exit(1)
    
    check_conda_env('torch')

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

    # ## Get all the subdirectoires of the base folder 
    # dir_list = next(os.walk(basecalled_folder))[1]
    # # print(intial_dir_list)

    given_dataset_list = sys.argv[2]
    if given_dataset_list == 'ALL_IN_DIR':
        print('All the folders inside the given folder will be considered.')
        dir_list = next(os.walk(basecalled_folder))[1]
    ## I can also make another option --exclude 
    else:
        dir_list = given_dataset_list.split(',')
        dir_list = list(filter(None, dir_list)) # remove empty strings 
        # print(dir_list)


    analysis_folder = sys.argv[3]
    analysis_folder = check_endWith(analysis_folder)

    thread = sys.argv[4]
    single_mult = sys.argv[5]
    len_thr = sys.argv[6]
    map_qual_thr = sys.argv[7]
    BP_FT_tougle = sys.argv[8]      ## it should have values: 'BP' or 'FT'



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


