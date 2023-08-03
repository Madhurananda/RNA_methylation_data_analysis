
import os,sys
sys.path.insert(0, '/home/madhu/work/codes')
from main import check_endWith, check_create_dir, check_env, write_file, count_lines

# import time, datetime
from datetime import datetime
import TS_finetune_MP

'''

This scripts does the model training. 

An example call: 
python /home/madhu/work/codes/ML_codes/madhu_codes/TS_finetune_run_MP.py /home/madhu/work/class_outputs/p5_decoding_epitranscriptional_native_RNA_seq /home/madhu/work/ref_transcriptome/IVT_seq/IVT_seq.fa model_output__train_p5_m6A_3_512_1500_512_1__A/ 1 3 512 6 F 1500 512 1 A A1 |& tee /home/madhu/Logs/finetune__train_p5_m6A_3_512_1500_512_1__A.log

Inputs are: 
basefolder, REF_GENOME, output_model_dir, TRAIN_epoch, layers, hidden_layers, cuda_id, F, lr, batch_size, downsampling_rate, METHYL_TYPE, TEST_CHR (optional)

'''

if __name__=='__main__' and '__file__' in globals():

    check_env('torch')


    startTime = datetime.now()
    current_time = startTime.strftime("%Y/%m/%d at %H:%M:%S")
    print('\n\nThe script starting at: ' + str(current_time), ' \n\n' )

    ###############################################################################
    ###### The following needs to be updated #######
    # basefolder = '/home/madhu/work/class_outputs/p3_m6A_RNA_modification_native_RNA_seq'
    # TRAIN_epoch = 1
    # REF_GENOME = '/mnt/labshare/share/reference_genome/EpiNano_Reference_sequences/cc.fasta'

    basefolder = sys.argv[1]
    REF_GENOME = sys.argv[2]
    output_model_dir = sys.argv[3]
    TRAIN_epoch = int(sys.argv[4])
    # METHYL_TYPE = sys.argv[12]

    ###############################################################################

    # basefolder = '/home/qliu10/projects/TAD/analysis/ft_tfrecordguppy_v5.0.11//';             #      revise
    # basefolder = '/home/madhu/work/analysis/p3_m6A_RNA_modification_native_RNA_seq'
    
    basefolder = check_endWith(basefolder)

    para_dict = {}
    para_dict['save_folder'] = basefolder + output_model_dir;  #      revise
    # para_dict['save_folder'] = '/home/qliu10/projects/TAD/analysis/testModelguppy_v5.0.11/';  #      revise

    check_create_dir(para_dict['save_folder'])
    # if os.path.isdir( para_dict['save_folder'] ):
    #     os.system('mkdir {}'.format( para_dict['save_folder'] ))

    para_dict['save_file_pref'] = 'BilstmMean.'+'layers'+sys.argv[5]+'.hs'+sys.argv[6]+'.'+sys.argv[8]+'.lr'+sys.argv[9]+'.b'+sys.argv[10]+(".p"+str(int(float(sys.argv[11])*1000+0.5)) if len(sys.argv)>7 else "")
    
    print(basefolder, para_dict['save_file_pref'])

    # 337,837
    para_dict['warmup_seq'] = 10000
    para_dict['save_frequency'] = 10000 ; #para_dict['warmup_seq']
    para_dict['log_frequency']  =   100 ; #para_dict['warmup_seq']//10
    
    # para_dict['train_epoch']  = 5
    para_dict['train_epoch']  = TRAIN_epoch;                                                          #      revise
    para_dict['size_layers'] =int(sys.argv[5])
    para_dict['size_hidden'] =int(sys.argv[6])
    para_dict['METHYL_TYPE'] = sys.argv[12]

    print('The methyl type is: ', para_dict['METHYL_TYPE'])

    # para_dict['cpu_devices'] = 3; #
    para_dict['cuda_devices'] = 1;
    para_dict['cuda_id'] = sys.argv[7];
    # print('para_dict: ', para_dict['cuda_id'])
    para_dict['adam_learning_rate'] = int(sys.argv[9])/1e7
    para_dict['input_batch_size'] = int(sys.argv[10])
    if len(sys.argv)>11:
        para_dict['downsampling_rate'] = float(sys.argv[11])
        if para_dict['downsampling_rate'] < 1e-10:
            del para_dict['downsampling_rate']
        if 'downsampling_rate' in para_dict and para_dict['downsampling_rate']<0.03:
            para_dict['train_epoch']  = 1000
        del para_dict['downsampling_rate']
    
    if len(sys.argv)>14:
        para_dict['ft_learning_rate'] = float(sys.argv[14])
        para_dict['save_file_pref'] = para_dict['save_file_pref'] + ".ft"+str(int(float(sys.argv[14])*1e7+0.5))
    if len(sys.argv)>15:
        para_dict['saved_model'] = sys.argv[15]

    # for DNA;
    #para_dict['length_thr'] = 1000
    # for RNA ## run on 2023/03/04 forget to comment this                                       #      revise
    para_dict['length_thr'] = 200

    # print(datetime.datetime.now())
    try:
        is_ft_t = (True if sys.argv[8] in [1, '1', 'T', 'True', 'true'] else False)
        ## The first one will have the label 0 and the second one 1. 
        #datafs = [basefolder+'/FAH58492_MN17479/', basefolder+'/MN17273_FAH58548/']
        #datafs = [basefolder+'/MN17273_FAH58548/', basefolder+'/FAH58492_MN17479/']

        # datafs = [basefolder+'/GSM3528749/TFRec_bp_FT/', basefolder+'GSM3528750/TFRec_bp_FT/']                 #      revise
        datafs = [basefolder+'train/class_0/', basefolder+'train/class_1/']                 #      revise
        
        if len(sys.argv)>13:
            if len(sys.argv[13]) == 2:
                TS_finetune_MP.finetuneTrainer( para_dict = para_dict, datafolders=datafs, is_ft=is_ft_t, ref_seq_file=REF_GENOME, index_file="tfrecord.index", random_seed=3, not_use_chr=sys.argv[13] )
            else:
                # print('(sys.argv[13]: ', sys.argv[13])
                # print('len(sys.argv[13]): ', len(sys.argv[13]))
                print('\n\n\n****** The chromosome not to be used should only have the length: 2, such as U1 *******\n\n\n')
        else:
            TS_finetune_MP.finetuneTrainer( para_dict = para_dict, datafolders=datafs, is_ft=is_ft_t, ref_seq_file=REF_GENOME, index_file="tfrecord.index", random_seed=3 );          #      revise
        
        # TS_finetune_MP.finetuneTrainer( para_dict = para_dict, datafolders=datafs, is_ft=is_ft_t, ref_seq_file=REF_GENOME, index_file="tfrecord.index", random_seed=3, not_use_chr=['cc6m_2244_T7_ecorv'] )


    except Exception as e:
        print(traceback.format_exc())
        print(e)
        # print(datetime.datetime.now())
    # print(datetime.datetime.now())

    executionTime = (datetime.now() - startTime)
    current_time = datetime.now().strftime("%Y/%m/%d at %H:%M:%S")
    print('\n\nThe script completed at: ' + str(current_time))
    print('Execution time: ' + str(executionTime), ' \n\n')


