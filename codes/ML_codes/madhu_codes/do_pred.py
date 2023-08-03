import os
import sys
sys.path.insert(0, '/home/madhu/work/codes')

from main import check_endWith, check_create_dir, check_env



'''
This script will do training and create log files etc. 

python /home/madhu/work/codes/ML_codes/madhu_codes/TS_finetune_pred_run_MP.py /home/madhu/work/class_outputs/p3_m6A_RNA_modification_native_RNA_seq/ /home/madhu/work/class_outputs/p3_m6A_RNA_modification_native_RNA_seq/model_PRED_output__train_p3_ALL_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.10258.pt__A/ 3 512 5 256 /home/madhu/work/class_outputs/p3_m6A_RNA_modification_native_RNA_seq/model_output__train_p3_3_512_1000_256_1__A/BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.10258.pt /mnt/labshare/share/reference_genome/EpiNano_Reference_sequences/cc.fasta /home/madhu/work/class_outputs/p3_m6A_RNA_modification_native_RNA_seq A |& tee /home/madhu/Logs/predict__train_p3_ALL_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.10258.pt__A.log

The inputs are: 
basefolder, basefolder+output_model_dir, n_layer, a_n_N_each_layer, gpu_ID, batch_size, basefolder+saved_model_dir+'/'+saved_model, REF_GENOME, pred_folder, METHYL_TYPE, logF_name

'''



if __name__=='__main__':
    

    # ## Test on p5 (IVT data)
    # REF_GENOME = '/home/madhu/work/ref_transcriptome/IVT_seq/IVT_seq.fa'
    # pred_folder = '/home/madhu/work/class_outputs/p5_decoding_epitranscriptional_native_RNA_seq'

    # ## train and test on p3 (seperate datasets)
    # REF_GENOME = '/mnt/labshare/share/reference_genome/EpiNano_Reference_sequences/cc.fasta'
    # pred_folder = '/home/madhu/work/class_outputs/p3_m6A_RNA_modification_native_RNA_seq'
    
    given_DIR = sys.argv[1]
    given_DIR = check_endWith(given_DIR)
    saved_model_dir = given_DIR.split('/')[len(given_DIR.split('/'))-2]
    basefolder = given_DIR.replace(saved_model_dir+'/', '')
    saved_model = sys.argv[2]
    gpu_ID = sys.argv[3]
    pred_token = sys.argv[4]
    REF_GENOME = sys.argv[5]
    pred_folder = sys.argv[6]

    basefolder = check_endWith(basefolder)

    METHYL_TYPE = saved_model_dir.split('__')[2]
    downsampling_rate = saved_model_dir.split('__')[1].split('_')[len(saved_model_dir.split('__')[1].split('_'))-1]
    batch_size = 1024   # I finally decided to run 1 job at a time with the max possible (1024) batch size 
    # batch_size = 256   # For running 3 multiple predictions parallelly
    # batch_size = saved_model_dir.split('__')[1].split('_')[len(saved_model_dir.split('__')[1].split('_'))-2]
    learning_rate = saved_model_dir.split('__')[1].split('_')[len(saved_model_dir.split('__')[1].split('_'))-3]
    a_n_N_each_layer = saved_model_dir.split('__')[1].split('_')[len(saved_model_dir.split('__')[1].split('_'))-4]
    n_layer = saved_model_dir.split('__')[1].split('_')[len(saved_model_dir.split('__')[1].split('_'))-5]


    output_model_dir = 'model_PRED_output__{}_{}_{}__{}/'.format(pred_token, saved_model_dir.split('__')[-1], saved_model, METHYL_TYPE)
    # print('output_model_dir: ', output_model_dir)

    logF_name = '/home/madhu/Logs/predict__{}_{}_{}__{}.log'.format(pred_token, saved_model_dir.split('__')[-1], saved_model, METHYL_TYPE) 

    predict_cmd = 'python /home/madhu/work/codes/ML_codes/madhu_codes/TS_finetune_pred_run_MP.py {} {} {} {} {} {} {} {} {} {} |& tee {}'.format(basefolder, basefolder+output_model_dir, n_layer, a_n_N_each_layer, gpu_ID, batch_size, basefolder+saved_model_dir+'/'+saved_model, REF_GENOME, pred_folder, METHYL_TYPE, logF_name)

    print( predict_cmd)
    print('\n*** ***\n')

