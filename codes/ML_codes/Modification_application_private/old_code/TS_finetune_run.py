
import os,sys
import time, datetime

import TS_finetune


if __name__=='__main__' and '__file__' in globals():
    para_dict = {}
    para_dict['save_folder'] = '/home/qliu10/projects/TAD/analysis/testModelguppy_v5.0.11/';  #      revise
    if os.path.isdir( para_dict['save_folder'] ):
        os.system('mkdir {}'.format( para_dict['save_folder'] ))
    para_dict['save_file_pref'] = 'BilstmMean.'+'layers'+sys.argv[1]+'.hs'+sys.argv[2]+'.'+sys.argv[4]+'.lr'+sys.argv[5]+'.b'+sys.argv[6]+(".p"+str(int(float(sys.argv[7])*1000+0.5)) if len(sys.argv)>7 else "")
    basefolder = '/home/qliu10/projects/TAD/analysis/ft_tfrecordguppy_v5.0.11//';             #      revise
    basefolder = '/home/madhu/datasets/analysis/p3_m6A_RNA_modification_native_RNA_seq/'
    print(basefolder, para_dict['save_file_pref'])

    # 337,837
    para_dict['warmup_seq'] = 10000
    para_dict['save_frequency'] = 10000 ; #para_dict['warmup_seq']
    para_dict['log_frequency']  =   100 ; #para_dict['warmup_seq']//10

    para_dict['train_epoch']  = 5
    para_dict['train_epoch']  = 1;                                                          #      revise
    para_dict['size_layers'] =int(sys.argv[1])
    para_dict['size_hidden'] =int(sys.argv[2])

    para_dict['cpu_devices'] = 3; #
    para_dict['cuda_devices'] = 1;
    para_dict['cuda_id'] = sys.argv[3];

    para_dict['adam_learning_rate'] = int(sys.argv[5])/1e7
    para_dict['input_batch_size'] = int(sys.argv[6])
    if len(sys.argv)>7:
        para_dict['downsampling_rate'] = float(sys.argv[7])
        if para_dict['downsampling_rate'] < 1e-10:
            del para_dict['downsampling_rate']
        if 'downsampling_rate' in para_dict and para_dict['downsampling_rate']<0.03:
            para_dict['train_epoch']  = 1000
        del para_dict['downsampling_rate']
    if len(sys.argv)>8:
        para_dict['ft_learning_rate'] = float(sys.argv[8])
        para_dict['save_file_pref'] = para_dict['save_file_pref'] + ".ft"+str(int(float(sys.argv[8])*1e7+0.5))
    if len(sys.argv)>9:
        para_dict['saved_model'] = sys.argv[9]

    # for DNA;
    #para_dict['length_thr'] = 1000
    # for RNA ## run on 2023/03/04 forget to comment this                                       #      revise
    para_dict['length_thr'] = 200

    print(datetime.datetime.now())
    try:
        is_ft_t = (True if sys.argv[4] in [1, '1', 'T', 'True', 'true'] else False)
        #datafs = [basefolder+'/FAH58492_MN17479/', basefolder+'/MN17273_FAH58548/']
        #datafs = [basefolder+'/MN17273_FAH58548/', basefolder+'/FAH58492_MN17479/']

        datafs = [basefolder+'/GSM3528751/TFRec_bp/', basefolder+'GSM3528752/TFRec_bp/:']                 #      revise

        #finetuneTrainer( para_dict = para_dict, datafolders=datafs, is_ft=is_ft_t, ref_seq_file='EpiNano/Reference_sequences/cc.fasta',index_file="tfrecord.index", adam_learning_rate=int(sys.argv[5])/1e7, input_batch_size=int(sys.argv[6]), random_seed=3, not_use_chr=['cc6m_2244_T7_ecorv'] , downsampling_rate=(float(sys.argv[7])if len(sys.argv)>7 else None), saved_model=(sys.argv[9] if len(sys.argv)>9 else None))
        TS_finetune.finetuneTrainer( para_dict = para_dict, datafolders=datafs, is_ft=is_ft_t, ref_seq_file='/mnt/labshare/share/reference_genome/yeast/sacCer3.fa',index_file="tfrecord.index", random_seed=3 );          #      revise
    except Exception as e:
        print(traceback.format_exc())
        print(e)
        print(datetime.datetime.now())
    print(datetime.datetime.now())


