
import os,sys
import time, datetime

import TS_finetune_pred




if __name__=='__main__' and '__file__' in globals():
    para_dict = {}
    # para_dict['save_folder'] = '/home/qliu10/projects/TAD/analysis/testPredguppy_v5.0.11/';
    para_dict['save_folder'] = '/home/madhu/work/class_outputs/p3_m6A_RNA_modification_native_RNA_seq/model_PRED_output/'
    if not os.path.isdir( para_dict['save_folder'] ):
        os.system( 'mkdir -p {}'.format( para_dict['save_folder'] ) )
    #para_dict['save_file_pref'] = 'BilstmMean.'+'layers'+sys.argv[1]+'.hs'+sys.argv[2]+'.'+sys.argv[4]+'.lr'+sys.argv[5]+'.b'+sys.argv[6]+'.'
    # basefolder = '/home/qliu10/projects/TAD/analysis/ft_tfrecordguppy_v5.0.11/';
    basefolder = '/home/madhu/work/analysis/p3_m6A_RNA_modification_native_RNA_seq/'

    saved_model = sys.argv[5]

    # analysis/testModelguppy_v5.0.11/BilstmMean.layers3.hs1024.F.lr500.b32.p1000.GPU.ep1.10025.pt
    #para_dict['save_file'] = 'analysis/pred_resultguppy_v5.0.11/'+'ps.'+'layers'+sys.argv[1]+'.hs'+sys.argv[2]+'.'+sys.argv[4]+'.predres'
    # para_dict['save_file'] = 'analysis/pred_resultguppy_v5.0.11/'+'ps.'+ '.'.join( saved_model.split('/')[-1].split('.')[1:-1] ) +'.predres' 
    para_dict['save_file'] = '/home/madhu/work/class_outputs/p3_m6A_RNA_modification_native_RNA_seq/model_PRED_output/'+'ps.'+ '.'.join( saved_model.split('/')[-1].split('.')[1:-1] ) +'.predres' 

    print(basefolder, para_dict['save_file'])

    para_dict['size_layers'] =int(sys.argv[1])
    para_dict['size_hidden'] =int(sys.argv[2])

    para_dict['cpu_devices'] = 15; #
    para_dict['cuda_devices'] = 1;
    para_dict['cuda_id'] = sys.argv[3];

    # for DNA;
    # para_dict['length_thr'] = 1000
    # for RNA
    para_dict['length_thr'] = 200

    
    print(datetime.datetime.now())
    try:
        ## The first one will have the label 0 and the second one 1. 
        
        #datafs = [basefolder+'/tfrecord_ft.5M/Curlcake_m6a_cc6m_2244_T7_ecorv/', basefolder+'/tfrecord_ft.5M/Curlcake_non_cc6m_2244_T7_ecorv/']
        # datafs = [basefolder+'/MN17273_FAH58548/', basefolder+'/FAH58492_MN17479/']
        # datafs = [basefolder+'/GSM3528749/TFRec_bp/', basefolder+'GSM3528750/TFRec_bp/:']
        datafs = [basefolder+'GSM3528750/TFRec_FT/', basefolder+'/GSM3528749/TFRec_FT/']
        TS_finetune_pred.pred( para_dict = para_dict, datafolders=datafs, saved_model=saved_model, ref_seq_file='/mnt/labshare/share/reference_genome/EpiNano_Reference_sequences/cc.fasta',index_file="tfrecord.index", input_batch_size=int(sys.argv[4]) )
    except Exception as e:
        print(traceback.format_exc())
        print(e)
        print(datetime.datetime.now())
    print(datetime.datetime.now())




