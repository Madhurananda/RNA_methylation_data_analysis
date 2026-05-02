
import os
import sys
sys.path.insert(0, '/home/madhu/work/codes')

from main import check_endWith, check_create_dir, check_env
from datetime import datetime
from glob import glob
from tqdm import tqdm
from multiprocessing.pool import ThreadPool, Pool
from multiprocessing import cpu_count
import numpy as np
import glob
from natsort import natsorted



if __name__=='__main__':
#     if len(sys.argv)<2:
#         print("Usage: {} {}".format(sys.argv[0], "basefolder"))
#         print("Example:")
#         print("       python {} {} {} {}".format(sys.argv[0], "/home/madhu/work/download/m6A_RNA_modification_native_RNA_seq/GSM3528751/extracted_data"))
# #           print("       python {} {} {} {}".format(sys.argv[0], "Curlcake_non"))
#         sys.exit(1)

    # check_conda_env('torch')

    '''
    I will execute something like this:
    python /home/madhu/work/codes/ML_codes/madhu_codes/do_pred_run.py /home/madhu/work/class_outputs/p3_m6A_RNA_modification_native_RNA_seq/model_output__train_p3_BOTH_3_512_1000_256_1__A /home/madhu/work/ref_transcriptome/IVT_seq/IVT_seq.fa /home/madhu/work/class_outputs/p5_decoding_epitranscriptional_native_RNA_seq X train_p3_BOTH_pred_p5_m6A
    '''


    model_skip = 5  # higher number means less number of models to be selected ... 



    saved_model_folder = sys.argv[1]
    REF_GENOME = sys.argv[2]
    pred_folder = sys.argv[3]
    # saved_model_folder = '/home/madhu/work/class_outputs/p3_m6A_RNA_modification_native_RNA_seq/model_output__train_p3_3_512_1000_256_1'
    list_GPU_IDs = sys.argv[4].split(',')
    pred_token = sys.argv[5]
    # N_parallel_jobs = int(sys.argv[3])
    N_parallel_jobs = 1

    saved_model_folder = check_endWith(saved_model_folder)
    
    # list_saved_model = glob.glob(saved_model_folder+'*.pt')
    list_saved_model = natsorted([os.path.join(root, name)
                    for root, dirs, files in os.walk(saved_model_folder)
                    for name in files
                    if name.endswith((".pt"))])
    
    list_saved_model = list_saved_model[0::model_skip]

    # print(len(list_saved_model[0::3]))
    # splited_models_list = list(np.array_split(list_saved_model, int(len(list_saved_model)/N_parallel_jobs)))
    # print('Len splited_models_list: ', len(splited_models_list))

    splited_models_list = []
    for i in range(0, len(list_saved_model), N_parallel_jobs):
        splited_models_list.append(list_saved_model[i:i+N_parallel_jobs])
    # print('Len splited_models_list: ', len(splited_models_list))


    i_GPU = 0
    print('\n\n')
    for models in splited_models_list:
        # print(models, 'on GPU: ', list_GPU_IDs[i_GPU])
        for a_model in models:
            # logF_name = '/home/madhu/Logs/ALL_doPREDrun_predict__{}_{}.log'.format(pred_token, a_model.split('/')[-1]) 
            pred_CMD = 'python /home/madhu/work/codes/ML_codes/madhu_codes/do_pred.py {} {} {} {} {} {}'.format(saved_model_folder, a_model.split('/')[-1], list_GPU_IDs[i_GPU], pred_token, REF_GENOME, pred_folder)
            # pred_CMD = 'python /home/madhu/work/codes/ML_codes/madhu_codes/do_pred.py {} {} {} {} > {} 2>&1'.format(saved_model_folder, a_model.split('/')[-1], list_GPU_IDs[i_GPU], pred_token, logF_name)
            # print( pred_CMD, '&\n')
            
            # print( pred_CMD)
            # print('\n*** ***\n')

            sys.stdout.flush()
            os.system(pred_CMD)
        # print('wait')

        i_GPU += 1
        if i_GPU == len(list_GPU_IDs):
            i_GPU = 0
    
    print('*** ***\n')


