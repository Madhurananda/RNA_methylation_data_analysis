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


Update later .... 

'''



if __name__=='__main__':
    if len(sys.argv)<2:
        print("Usage: {} {}".format(sys.argv[0], "basefolder"))
        print("Example:")
        print("       python {} {} {} {}".format(sys.argv[0], "/home/madhu/work/download/m6A_RNA_modification_native_RNA_seq/GSM3528751/extracted_data"))
#           print("       python {} {} {} {}".format(sys.argv[0], "Curlcake_non"))
        sys.exit(1)
    
    # check_env('torch')

    # startTime = datetime.now()
    # current_time = startTime.strftime("%Y/%m/%d at %H:%M:%S")
    # print('\n\nThe script starting at: ' + str(current_time), ' \n\n' )

    '''
    Example call: 
    python /home/madhu/work/codes/ML_codes/madhu_codes/do_train.py train_p5_m6A A1
    python /home/madhu/work/codes/ML_codes/madhu_codes/do_train.py train_p5_m5C C1


    I will execute something like this, generated from this script. 
    python /home/madhu/work/codes/ML_codes/madhu_codes/TS_finetune_run_MP.py 3 512 1 F 500 32 1 |& tee /home/madhu/Logs/finetune_TEST_l2.log 
    '''

    training_token =sys.argv[1]



    ######################## Training on paper 3 #########################
    # basefolder = '/home/madhu/work/class_outputs/p3_m6A_RNA_modification_native_RNA_seq'
    # REF_GENOME = '/mnt/labshare/share/reference_genome/EpiNano_Reference_sequences/cc.fasta'
    # METHYL_TYPE = 'A'  ## Can have values: A, C, G, U
    ###########################################################################


    ######################## Training on paper 3 and paper 5 #########################
    # basefolder = '/home/madhu/work/class_outputs/p5_decoding_epitranscriptional_native_RNA_seq'
    # REF_GENOME = '/home/madhu/work/ref_transcriptome/cc_IVT/cc_IVT.fa'
    # METHYL_TYPE = 'A'  ## Can have values: A, C, G, U
    ###########################################################################



    ######################## Training on paper 5 (IVT data) #########################
    '''
    A:
    m1A
    m6A
    IVT_normalA_fast5

    C:
    hm5C
    FormylC
    m5C
    IVT_normalC_fast5

    G:
    m7G
    Inosine
    IVT_normalG_fast5

    U:
    pseudoU
    IVT_normalU_fast5

    chr_list = ['A1', 'A2', 'U1', 'U2', 'C1', 'C2', 'G1', 'G2']
    '''
    basefolder = '/home/madhu/work/class_outputs/p5_decoding_epitranscriptional_native_RNA_seq'
    REF_GENOME = '/home/madhu/work/ref_transcriptome/IVT_seq/IVT_seq.fa'
    ###########################################################################
    


    ###########################################################################
    # ## Training: pseudoU (U) 
    # METHYL_TYPE = 'U'  ## Can have values: A, C, G, U 
    # # TEST_NOT_use_chr = ['U1']
    # if len(sys.argv) >2:
    #     TEST_NOT_use_chr = sys.argv[2]
    
    ## Training: m6A (A)
    METHYL_TYPE = 'A'  ## Can have values: A, C, G, U 
    if len(sys.argv) >2:
        TEST_NOT_use_chr = sys.argv[2]
    
    # ## Training: m5C (C)
    # METHYL_TYPE = 'C'  ## Can have values: A, C, G, U 
    # if len(sys.argv) >2:
    #     TEST_NOT_use_chr = sys.argv[2]
    ###########################################################################



    ###########################################################################
    ## Set up the parameters ... 
    n_layers = [3]
    n_N_each_layer = [1024]
    learning_rates = [1500]
    batch_sizes = [1024]
    TRAIN_epoch = 3
    gpu_ID = 'X'
    ###########################################################################



    downsampling_rates = [1]

    print('\n\n')
    for n_layer in n_layers:
        for a_n_N_each_layer in n_N_each_layer:
            for learning_rate in learning_rates:
                for batch_size in batch_sizes:
                    for downsampling_rate in downsampling_rates:
                        # train_cmd = 'python /home/madhu/work/codes/ML_codes/madhu_codes/TS_finetune_run_MP.py 3 512 1 F 500 32 1 |& tee /home/madhu/Logs/finetune_TEST_l2.log'
                        # NEW_train_cmd = 'python /home/madhu/work/codes/ML_codes/madhu_codes/TS_finetune_run_MP.py /home/madhu/work/class_outputs/p3_m6A_RNA_modification_native_RNA_seq /mnt/labshare/share/reference_genome/EpiNano_Reference_sequences/cc.fasta output_model_dir TRAIN_epoch 3 512 1 F 500 32 1 |& tee /home/madhu/Logs/finetune_TEST_l2.log'

                        output_model_dir = 'model_output__{}_{}_{}_{}_{}_{}_{}__{}/'.format(training_token, n_layer, a_n_N_each_layer,learning_rate, batch_size, downsampling_rate, TRAIN_epoch, METHYL_TYPE)

                        # print('n_layer: ', n_layer)
                        # print('a_n_N_each_layer: ', a_n_N_each_layer)
                        # print('GPU ID: ', gpu_ID)
                        # print('learning_rate: ', learning_rate)
                        # print('batch_size: ', batch_size)
                        # print('downsampling_rate: ', downsampling_rate)
                        # print('output_model_dir: ', output_model_dir)

                        logF_name = '/home/madhu/Logs/finetune__{}_{}_{}_{}_{}_{}_{}__{}.log'.format(training_token, n_layer, a_n_N_each_layer, learning_rate, batch_size, downsampling_rate, TRAIN_epoch, METHYL_TYPE) 
                        # train_cmd = 'python /home/madhu/work/codes/ML_codes/madhu_codes/TS_finetune_run_MP.py {} {} {} {} {} {} {} F {} {} {}'.format(basefolder, REF_GENOME, output_model_dir, TRAIN_epoch, n_layer, a_n_N_each_layer, gpu_ID, learning_rate, batch_size, downsampling_rate)

                        if 'TEST_NOT_use_chr' in locals():
                            train_cmd = 'python /home/madhu/work/codes/ML_codes/madhu_codes/TS_finetune_run_MP.py {} {} {} {} {} {} {} F {} {} {} {} {} |& tee {}'.format(basefolder, REF_GENOME, output_model_dir, TRAIN_epoch, n_layer, a_n_N_each_layer, gpu_ID, learning_rate, batch_size, downsampling_rate, METHYL_TYPE, TEST_NOT_use_chr, logF_name)
                        else:
                            train_cmd = 'python /home/madhu/work/codes/ML_codes/madhu_codes/TS_finetune_run_MP.py {} {} {} {} {} {} {} F {} {} {} {} |& tee {}'.format(basefolder, REF_GENOME, output_model_dir, TRAIN_epoch, n_layer, a_n_N_each_layer, gpu_ID, learning_rate, batch_size, downsampling_rate, METHYL_TYPE, logF_name)

                        # train_cmd = 'python /home/madhu/work/codes/ML_codes/madhu_codes/TS_finetune_run_MP.py {} {} {} {} {} {} {} F {} {} {} |& tee {}'.format(basefolder, REF_GENOME, output_model_dir, TRAIN_epoch, n_layer, a_n_N_each_layer, gpu_ID, learning_rate, batch_size, downsampling_rate, logF_name)
                        
                        # print('logF_name: ', logF_name)
                        # print('train_cmd is: ')
                        print(train_cmd)

                        # sys.stdout.flush()
                        # os.system(train_cmd)
                        # stream = os.popen(train_cmd)
                        # output = stream.read()
                        # print('\n\nThe Training output: ------------------------------------------------ \n', (output), '\n\n')

                    print('\n*** ***\n')


    # executionTime = (datetime.now() - startTime)
    # current_time = datetime.now().strftime("%Y/%m/%d at %H:%M:%S")
    # print('\n\nThe script completed at: ' + str(current_time))
    # print('Execution time: ' + str(executionTime), ' \n\n')


