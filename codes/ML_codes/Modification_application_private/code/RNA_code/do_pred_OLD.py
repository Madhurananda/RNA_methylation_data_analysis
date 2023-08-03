import os
import sys
sys.path.insert(0, '/home/madhu/work/codes')

from main import check_endWith, check_create_dir, check_env
from datetime import datetime
from glob import glob
from tqdm import tqdm
from multiprocessing.pool import ThreadPool, Pool
from multiprocessing import cpu_count


'''
This script will do training and create log files etc. 



'''



if __name__=='__main__':
#     if len(sys.argv)<2:
#         print("Usage: {} {}".format(sys.argv[0], "basefolder"))
#         print("Example:")
#         print("       python {} {} {} {}".format(sys.argv[0], "/home/madhu/work/download/m6A_RNA_modification_native_RNA_seq/GSM3528751/extracted_data"))
# #           print("       python {} {} {} {}".format(sys.argv[0], "Curlcake_non"))
#         sys.exit(1)

    # check_env('torch')

    # startTime = datetime.now()
    # current_time = startTime.strftime("%Y/%m/%d at %H:%M:%S")
    # print('\n\nThe script starting at: ' + str(current_time), ' \n\n' )

    '''
    I will execute something like this:
    python /home/madhu/work/codes/ML_codes/madhu_codes/TS_finetune_pred_run_MP.py 3 512 0 32 /home/madhu/work/analysis/p3_m6A_RNA_modification_native_RNA_seq/model_output/BilstmMean.layers3.hs512.F.lr500.b32.p1000.GPU.ep1.40543.pt |& tee /home/madhu/LOGs/finetune_pred_TEST_l2.log 
    '''

    # basefolder = '/home/madhu/work/class_outputs/p3_m6A_RNA_modification_native_RNA_seq'
    

    # ## train on p3, test on p5
    # REF_GENOME = '/home/madhu/work/ref_transcriptome/IVT_seq/IVT_seq.fa'
    # pred_folder = '/home/madhu/work/class_outputs/p5_decoding_epitranscriptional_native_RNA_seq'

    ## train and test on p3 (seperate datasets)
    REF_GENOME = '/mnt/labshare/share/reference_genome/EpiNano_Reference_sequences/cc.fasta'
    pred_folder = '/home/madhu/work/class_outputs/p3_m6A_RNA_modification_native_RNA_seq'
    
    # saved_model_dir = 'model_output__train_p3_3_512_1000_256_1'
    # saved_model_dir = 'model_output__train_p3_3_1024_1500_512_1__A'
    # saved_model = 'BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.550227.pt'
    # gpu_ID = 0

    given_DIR = sys.argv[1]
    given_DIR = check_endWith(given_DIR)
    saved_model_dir = given_DIR.split('/')[len(given_DIR.split('/'))-2]
    basefolder = given_DIR.replace(saved_model_dir+'/', '')

    # saved_model_dir = sys.argv[1]
    saved_model = sys.argv[2]
    gpu_ID = sys.argv[3]
    pred_token = sys.argv[4]

    basefolder = check_endWith(basefolder)

    # downsampling_rate = saved_model_dir.split('_')[len(saved_model_dir.split('_'))-1]
    # # batch_size = 1024  # The maximum number of batch size which can be loaded to a GPU
    # batch_size = 256   # For running multiple predictions parallelly 
    # # batch_size = saved_model_dir.split('_')[len(saved_model_dir.split('_'))-2]
    # learning_rate = saved_model_dir.split('_')[len(saved_model_dir.split('_'))-3]
    # a_n_N_each_layer = saved_model_dir.split('_')[len(saved_model_dir.split('_'))-4]
    # n_layer = saved_model_dir.split('_')[len(saved_model_dir.split('_'))-5]

    METHYL_TYPE = saved_model_dir.split('__')[2]
    downsampling_rate = saved_model_dir.split('__')[1].split('_')[len(saved_model_dir.split('__')[1].split('_'))-1]
    batch_size = 256   # For running 3 multiple predictions parallelly 
    # batch_size = saved_model_dir.split('__')[1].split('_')[len(saved_model_dir.split('__')[1].split('_'))-2]
    learning_rate = saved_model_dir.split('__')[1].split('_')[len(saved_model_dir.split('__')[1].split('_'))-3]
    a_n_N_each_layer = saved_model_dir.split('__')[1].split('_')[len(saved_model_dir.split('__')[1].split('_'))-4]
    n_layer = saved_model_dir.split('__')[1].split('_')[len(saved_model_dir.split('__')[1].split('_'))-5]



    # ## Find the saved models inside the folder ... 
    # model_lists = os.listdir(basefolder+saved_model_dir)
    # print('file_lists: ', model_lists)


    output_model_dir = 'model_PRED_output__{}_{}_{}__{}/'.format(pred_token, saved_model_dir.split('__')[-1], saved_model, METHYL_TYPE)
    # print('output_model_dir: ', output_model_dir)

    logF_name = '/home/madhu/LOGs/predict__{}_{}_{}__{}.log'.format(pred_token, saved_model_dir.split('__')[-1], saved_model, METHYL_TYPE) 

    predict_cmd = 'python /home/madhu/work/codes/ML_codes/madhu_codes/TS_finetune_pred_run_MP.py {} {} {} {} {} {} {} {} {} {} |& tee {}'.format(basefolder, basefolder+output_model_dir, n_layer, a_n_N_each_layer, gpu_ID, batch_size, basefolder+saved_model_dir+'/'+saved_model, REF_GENOME, pred_folder, METHYL_TYPE, logF_name)

    # predict_cmd = 'python /home/madhu/work/codes/ML_codes/madhu_codes/TS_finetune_pred_run_MP.py {} {} {} {} {} {} {} {} {} {} > {} 2>&1'.format(basefolder, basefolder+output_model_dir, n_layer, a_n_N_each_layer, gpu_ID, batch_size, basefolder+saved_model_dir+'/'+saved_model, REF_GENOME, pred_folder, METHYL_TYPE, logF_name)

    # print('logF_name: ', logF_name)
    # print('\n*********\npredict_cmd is: ', predict_cmd, '\n*********\n')
    print( predict_cmd)
    print('\n*** ***\n')
    # sys.stdout.flush()
    # os.system(predict_cmd)
    # stream = os.popen(predict_cmd)
    # output = stream.read()
    # print('\n\nThe Training output: ------------------------------------------------ \n', (output), '\n\n')


    # executionTime = (datetime.now() - startTime)
    # current_time = datetime.now().strftime("%Y/%m/%d at %H:%M:%S")
    # print('\n\nThe script completed at: ' + str(current_time))
    # print('Execution time: ' + str(executionTime), ' \n\n')


