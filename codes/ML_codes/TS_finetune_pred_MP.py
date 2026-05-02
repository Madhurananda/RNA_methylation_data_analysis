
import os,sys
import random
import time
import numpy as np
import copy
import math
import traceback

sys.path.insert(0, '/home/madhu/work/codes')
from main import save_file, load_file

import torch
import torch.multiprocessing

# from tqdm import tqdm

sys.path.append( os.path.dirname(__file__) )

#import TS_model
import TS_biLSTM
#import TS_dataset_reader

import TS_global
import TS_ftdataset_reader
import read_tombo

from datetime import datetime
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib.pyplot as plt



"""
best fine-tuning learning rate (among 5e-5, 4e-5, 3e-5, and 2e-5)

https://arxiv.org/pdf/2107.12460.pdf
Frozen pretrained: learning rate 1e-4, 1e-3;
Unfrozen pretrained: learning rate 1e-5, 1e-4

unfrozen is better than frozen
 After carefully redesigning the
empirical setup, we find that when tuning learning rates properly, pretrained transformers do outperform or match training from scratch in all of
our tasks, but only as long as the entire model is
fine-tuned. Thus, while transfer from pre-trained
language models to other modalities does indeed
provide gains and hints at exciting possibilities
for future work, properly tuning hyperparameters
is important for arriving at robust findings.

the best results are obtained when finetuning all of the weights of a pretrained model

 These results directly contradict prior work which
concluded that freezing most of the model leads to superior
performance. We demonstrate that these prior conclusions
were an artefact of using a specific, fixed learning rate
"""


class train_feat_Batch:
    def __init__(self, data):
        transposed_data = list(zip(*data))
        self.feat = torch.stack(transposed_data[0], 0)
        self.pad = torch.stack(transposed_data[1], 0)
        #self.modLab = torch.stack(transposed_data[2], 0)

    def pin_memory(self):
        self.feat = self.feat.pin_memory()
        self.pad = self.pad.pin_memory()
        #self.modLab = self.modLab.pin_memory()
        return self

def train_feat_Batch_wrapper(batch):
    return train_feat_Batch(batch)


def check_nucl(METHYL_TYPE):
    # print('The methyl type is: ', METHYL_TYPE)
    inITstr = ''
    N_inITstr = ''
    if METHYL_TYPE == 'A':
        # print('The methyl type is: ', METHYL_TYPE)
        inITstr = ['A', 'a']
        N_inITstr = ['T', 't', 'U', 'u']
    elif METHYL_TYPE == 'C':
        # print('The methyl type is: ', METHYL_TYPE)
        inITstr = ['C', 'c']
        N_inITstr = ['G', 'g']
    elif METHYL_TYPE == 'G':
        # print('The methyl type is: ', METHYL_TYPE)
        inITstr = ['G', 'g']
        N_inITstr = ['C', 'c']
    elif METHYL_TYPE == 'U':
        # print('The methyl type is: ', METHYL_TYPE)
        inITstr = ['T', 't', 'U', 'u']
        N_inITstr = ['A', 'a']
    else:
        print('PLease provide methyl types as: A, C, G, U')
    return inITstr, N_inITstr 


def generate_pred_feat(signal_full_info2, ref_seq_dict, seqr, para_dict):
    feat_ts_dict = {}
    pad_inf_dict = {}

    map_info = {}
    fn_list = []

    signal_full_length = copy.deepcopy(signal_full_info2[0])
    for _i_seq in range(len(signal_full_length)):
        #print(_i_seq, signal_full_length[_i_seq].size(1), para_dict['length_thr'])
        #print(signal_full_info2[3][_i_seq] )
        #print(sorted(list(ref_seq_dict.keys())))
        #if signal_full_length[_i_seq].size(1)<TS_global.min_seq_len: continue;
        if signal_full_length[_i_seq].size(1)<para_dict['length_thr']: continue;
        map_chr = signal_full_info2[3][_i_seq][0][0]
        if map_chr not in ref_seq_dict:
            continue;
        map_stn = signal_full_info2[3][_i_seq][1][0]
        map_sta = signal_full_info2[3][_i_seq][2].item()
        map_end = signal_full_info2[3][_i_seq][3].item()

        inITstr, N_inITstr = check_nucl(para_dict['METHYL_TYPE'])
        if inITstr == '' or N_inITstr == '':
            sys.exit("There are something wrong with METHYL TYPE.")
        
        a_pos_list = []
        r_pos_list = []
        for posi in range(map_end - map_sta):
            if map_stn=='+':
                #if ref_seq_dict[ map_chr ][ posi + map_sta ] in ['T', 't']: #['A', 'a']:
                if ref_seq_dict[ map_chr ][ posi + map_sta ] in inITstr:
                    a_pos_list.append( posi )
                    r_pos_list.append( posi + map_sta )
                    #print("test: {} {} {}/{} {} {} ".format( map_end, map_sta, map_stn,  posi, map_chr, posi + map_sta ))
            else:
                #if ref_seq_dict[ map_chr ][ map_end - posi - 1] in ['A', 'a']: # ['T', 't']:
                if ref_seq_dict[ map_chr ][ map_end - posi - 1] in N_inITstr:
                    a_pos_list.append( posi )
                    r_pos_list.append( map_end - posi - 1 )
                    #print("test: {} {} {}/{} {} {} ".format( map_end, map_sta, map_stn,  posi, map_chr, map_end - posi - 1 ) )

        use_ref_pos = []
        feat_ts = []
        pad_inf = []
        for _this_cp_ind in range(len(a_pos_list)):
            _this_cp = a_pos_list[_this_cp_ind]
            if _this_cp<5 or _this_cp>signal_full_length[_i_seq].size(1)-5: continue;
            use_ref_pos.append( r_pos_list[_this_cp_ind] )

            this_feat = []
            this_pad = []
            _this_start = _this_cp - (TS_global.input_len//2)
            _this_end = _this_cp + (TS_global.input_len//2)+1

            while _this_start<0:
                this_feat.append( torch.zeros( TS_global.dis_size ) )
                this_pad.append(0)
                _this_start +=  1;
            while _this_start < _this_end:
                if _this_start>=signal_full_length[_i_seq].size(1):
                    this_feat.append( torch.zeros( TS_global.dis_size ) )
                    this_pad.append(0)
                else:
                    this_feat.append( signal_full_length[_i_seq][0, _this_start, :] )
                    this_pad.append(1)
                _this_start +=  1;

            feat_ts.append( torch.stack(this_feat) )
            pad_inf.append( this_pad )

        if len( use_ref_pos )==0: continue;

        feat_ts_dict [ signal_full_info2[4][_i_seq][0] ] = feat_ts
        pad_inf_dict [ signal_full_info2[4][_i_seq][0] ] = pad_inf
        if signal_full_info2[4][_i_seq][0] in map_info:
            print("Duplicate files {} for {}: {}:{}:{}-{}".format( map_info[ signal_full_info2[4][_i_seq][0]],  signal_full_info2[4][_i_seq][0], map_chr,map_stn, map_sta, map_end  ))
        map_info[ signal_full_info2[4][_i_seq][0] ] = (map_chr,map_stn, map_sta, map_end, use_ref_pos);
        fn_list.append( signal_full_info2[4][_i_seq][0]  )
        
    return (feat_ts_dict, pad_inf_dict, map_info, fn_list );

def calculate_perf( cm ):
    num_true = cm[0][0] + cm[0][1];
    num_fals = cm[1][0] + cm[1][1];

    total = num_true + num_fals
    pred_t = cm[0][0] + cm[1][0]
    pred_f = cm[0][1] + cm[1][1]

    recall = cm[0][0]/(1 if num_true==0 else num_true)
    precision = cm[0][0]/(1 if pred_t==0 else pred_t)
    acc = (cm[0][0] + cm[1][1])/(1 if total==0 else total)
    f1 = 2*recall*precision/(recall + precision if recall + precision>0 else 1)

    print("R={:3f} P=={:3f} Acc={:3f} F1={:3f}".format( recall, precision, acc, f1 ))



def pred(para_dict, datafolders, ref_seq_file, saved_model, use_chr=[], index_file='tfrecord.index', input_batch_size=1024, cuda_devices: int=1, cpu_devices: int=5):
    torch.multiprocessing.set_sharing_strategy('file_system')

    if not os.path.isdir(para_dict['save_folder']):
        os.system('mkdir -p {}'.format(para_dict['save_folder']))

    if 'saved_model' in para_dict:
        saved_model = para_dict['saved_model']
    if 'cpu_devices' in para_dict:   
        if para_dict['cpu_devices']>0: cpu_devices   = para_dict['cpu_devices']
        else: print('cpu_devicesis too small: {}. not used.'.format(para_dict['cpu_devices']))
    if 'cuda_devices' in para_dict:  
        if para_dict['cuda_devices']>0: cuda_devices  = para_dict['cuda_devices']
        else:  print('cuda_devices too small: {}. not used.'.format(para_dict['cuda_devices']))
    if 'input_batch_size' in para_dict: 
        if para_dict['input_batch_size']>0: input_batch_size = para_dict['input_batch_size']
        else: 
            print('input_batch_sizeis too small: {}. not used.'.format(para_dict['input_batch_size']))
    if 'index_file' in para_dict:    
        index_file    = para_dict['index_file']

    print("Key learn info input_batch_size={} cudaID={} size_layers={} size_hidden={}".format(input_batch_size, ("cuda:"+para_dict['cuda_id'] if 'cuda_id' in para_dict else "cuda:0"), (3 if 'size_layers' not in para_dict or para_dict['size_layers']<1 else para_dict['size_layers']), (1024 if 'size_hidden' not in para_dict or para_dict['size_hidden']<10 else para_dict['size_hidden']) ))

    mod_signal_model = TS_biLSTM.SignalPred( ispretrain = False, pretrainmean=False,
                                                 size_layers = (3 if 'size_layers' not in para_dict or para_dict['size_layers']<1 else para_dict['size_layers']),
                                                 size_hidden = (1024 if 'size_hidden' not in para_dict or para_dict['size_hidden']<10 else para_dict['size_hidden'] ) );
    device = torch.device(("cuda:"+para_dict['cuda_id'] if 'cuda_id' in para_dict else "cuda:0") if torch.cuda.is_available() else "CPU")
    mod_signal_model = mod_signal_model.to(device);
    
    cpuCount = os.cpu_count()
    used_cpu_devices = cpu_devices
    if cpuCount <= used_cpu_devices : used_cpu_devices = cpuCount
    print("Number of CPUs in the system: {}. Used CPUs: {}".format(cpuCount, used_cpu_devices ) )

    seq_dict = TS_ftdataset_reader.read_seq_id( ref_seq_file )
    # print("All chrs={}".format( sorted(list(seq_dict.keys())) ) );
    # print("Use chrs={}".format( use_chr ) );
    if (not use_chr==None) and len(use_chr)>0:
        for nous_c in sorted(list(seq_dict.keys())):
            if nous_c not in use_chr:
                del seq_dict[ nous_c ]
    seqr = read_tombo.SeqReverse()

    if True:
        t_checkpoint = torch.load( saved_model ); #, map_location=device );
        mod_signal_model.load_state_dict(t_checkpoint['model_state_dict'])
        #_optmizer.load_state_dict(t_checkpoint['optimizer_state_dict'])
        num_train_seq = t_checkpoint['num_train_seq']
        save_epoch_i = t_checkpoint['epoch']
        train_epoch_i = t_checkpoint['epoch']-1
        loss_func = t_checkpoint['loss_func']
        num_train_1epoch = t_checkpoint['num_train_1epoch']
        max_estimat_seq = t_checkpoint['max_estimat_seq']
        total_max_estimat_seq = max_estimat_seq

    loss_func = torch.nn.CrossEntropyLoss()
    used_cpu_devices = 20

    start_time = time.time()
    pred_results = {}
    mod_signal_model.eval();
    loss_total = [0, 0]
    pred_size = [0, 0]

    read_pred_writer = open(para_dict['save_file'], 'w')

    ### Variables needed for ROC-AUC
    y_true_list = []
    y_score_list = []
    saveFig_READ_ROC_title = para_dict['save_folder']+para_dict['save_folder'].split('/')[len(para_dict['save_folder'].split('/'))-2]+'__READ_ROC_curve.png'


    process_i = 0;
    with torch.no_grad(): 
        for dataf_ind in range(len(datafolders)):
            dataf = datafolders[ dataf_ind ]
            fs_dataset_0 = TS_ftdataset_reader.TS_Dataset(dataf, index_file)
            fs_dataloader_0 = torch.utils.data.DataLoader(fs_dataset_0, batch_size=1, num_workers=used_cpu_devices, shuffle=True)
            # print(type(fs_dataset_0), len(fs_dataset_0), fs_dataset_0.total_file_seq_count, fs_dataset_0.total_seq_count )
            sys.stdout.flush()

            # print('\n\n******** fs_dataloader_0: ', len(fs_dataloader_0), '*********\n\n')
            testi = 0;
            batch_data_ind = 0
            for batch_fidx_0, batch_data_0 in enumerate( fs_dataloader_0 ):

                # gen_pred_feat_save_name = para_dict['save_folder']+'/'+str(dataf_ind)+'_'+str(batch_data_ind)+'__gen_pred_feat.pkl'
                # if not os.path.isfile(gen_pred_feat_save_name):
                #     feat_ts_dict, pad_inf_dict, map_info_dict, fn_list = generate_pred_feat ( batch_data_0, seq_dict, seqr, para_dict)
                #     save_file(gen_pred_feat_save_name, (feat_ts_dict, pad_inf_dict, map_info_dict, fn_list))
                #     print('generate_pred_feat Features were saved at: ', gen_pred_feat_save_name)
                # else:
                #     (feat_ts_dict, pad_inf_dict, map_info_dict, fn_list) = load_file(gen_pred_feat_save_name)
                #     print('generate_pred_feat was loaded from: ', gen_pred_feat_save_name)
                
                # print('\n\n******** batch_data_0: ', (batch_data_0), '*********\n\n')
                # print('\n\n******** seq_dict: ', (seq_dict), '*********\n\n')

                feat_ts_dict, pad_inf_dict, map_info_dict, fn_list = generate_pred_feat ( batch_data_0, seq_dict, seqr, para_dict)

                # print('feat_ts_dict: ', feat_ts_dict)
                batch_data_ind +=1

                # current_time = datetime.now().strftime("%Y/%m/%d at %H:%M:%S")
                # print('\n\ngenerate_pred_feat completed at: ' + str(current_time))
                
                ## fn_list is the list of fast5 reads 
                ## map_info[0], map_info[1] are: gene, strand (+ or -)
                ## map_info[4][ ref_pos_ind ] is: position 


                if len( fn_list ) ==0: continue;
                
                testi += 1; 
                #if testi >= fs_dataset_0.total_file_seq_count*0.1: break;

                for _fn_ in fn_list:
                    # if not _fn_ == 'GXB01170_20180726_FAH85615_GA10000_sequencing_run_RNAAB089716_65817_read_6637_ch_20_strand.fast5':
                    #     continue
                    # print('_fn_: ', _fn_, ' dataf_ind: ', dataf_ind)
                    process_i += 1;
                    map_info = map_info_dict[_fn_]

                    gene_id = map_info[0]
                    strand = map_info[1]

                    if gene_id not in pred_results:
                        pred_results[ gene_id ] = {};
                    if strand not in pred_results[ gene_id ]:
                        pred_results[ gene_id ][ strand ] = {};
                        
                    read_pred_writer.write(">{} {} {} {}\n".format( dataf_ind, _fn_, gene_id, strand ))

                    # print('len(feat_ts_dict[_fn_]): ', len(feat_ts_dict[_fn_]))

                    ###
                    for start_ind in range(0, len(feat_ts_dict[_fn_]), input_batch_size): 
                        end_ind = len(feat_ts_dict[_fn_]) if start_ind+input_batch_size*1.5>=len(feat_ts_dict[_fn_]) else start_ind+input_batch_size

                        # print('******* start_ind: ', start_ind)
                        # print('******* end_ind: ', end_ind)

                        t_feat = torch.stack(feat_ts_dict[_fn_][start_ind:end_ind], 0)
                        
                        # print('t_feat: ', t_feat)
                        # print('t_feat.shape', t_feat.shape)

                        # ## See, if there is any non-zero element in the tensor 
                        # print('len(torch.nonzero( t_feat )): ', len(torch.nonzero( t_feat )))
                        # if len(torch.nonzero( t_feat )) != 0:
                        #     print('t_feat: ', t_feat)
                        #     print('t_feat.shape', t_feat.shape)
                        # else:
                        #     print('FEATURES HAVE ONLY ZEROS .... ')


                        pred_mask = mod_signal_model.forward( t_feat.to(device), device );
                        if dataf_ind==0:
                            class_lab = torch.zeros(t_feat.size(0), dtype=torch.long ).to(device)
                        else: 
                            class_lab = torch.ones( t_feat.size(0), dtype=torch.long ).to(device)
                        
                        if len(pred_mask.size()) == 1:
                            print('\n\n******* pred_mask.size: ', pred_mask.size(), '********* \n\n')
                        
                        # ## This is for the site-based performance ... 
                        # y_true = np.array(class_lab.cpu())
                        # y_score = np.array(pred_mask.cpu())[:, 1]
                        # y_true_list.extend(y_true)
                        # y_score_list.extend(y_score)

                        loss = loss_func(pred_mask, class_lab)
                        loss_total[ dataf_ind ] += loss*t_feat.size(0)
                        pred_size[ dataf_ind ] += t_feat.size(0)

                        pred_id = torch.max(pred_mask, 1)[1]
                        pred_proba = pred_mask[:, 1]

                        # print('\n\n\n************** pred_mask shape: ', pred_proba.shape)
                        # print('pred_id shape: ', pred_id.shape, ' *************\n\n')

                        for ref_pos_ind in range(start_ind, end_ind):    
                            ref_pos = map_info[4][ ref_pos_ind ]
                            pred_label = pred_id[ref_pos_ind-start_ind].item()
                            pred_scores = pred_proba[ref_pos_ind-start_ind].item()

                            # ## This is for the read-based performance ... 
                            # y_true = [dataf_ind]
                            # y_score = [pred_scores]
                            # # print('\n\n***y_true: ', y_true)
                            # # print('y_score: ', y_score, '****\n\n')
                            # y_true_list.extend(y_true)
                            # y_score_list.extend(y_score)

                            read_pred_writer.write("{} {}\n".format( ref_pos, pred_label ))
                            if ref_pos not in pred_results[ gene_id ][ strand ]:
                                pred_results[ gene_id ][ strand ][ ref_pos ] = [[0, 0], [0, 0]]
                            pred_results[ gene_id ][ strand ][ ref_pos ][ dataf_ind ] [ pred_label ] += 1
                        if end_ind==len(feat_ts_dict[_fn_]):
                            break;
                    ###
                    read_pred_writer.write(">\n")

                    if process_i%5000==0:
                        print("{} Processing {} for {:.1f} {:.5f}/{:,} {}".format(dataf_ind, process_i, time.time()-start_time, loss_total[ dataf_ind ]/(pred_size[ dataf_ind ] if pred_size[ dataf_ind ]>0 else 1), pred_size[ dataf_ind ], t_feat.shape))
                        start_time = time.time()
                        sys.stdout.flush()
            
            print("{} Processing {} for {:.1f} {:.5f}/{:,} ".format(dataf_ind, process_i, time.time()-start_time, loss_total[ dataf_ind ]/(pred_size[ dataf_ind ] if pred_size[ dataf_ind ]>0 else 1), pred_size[ dataf_ind ] ))
            start_time = time.time()
            sys.stdout.flush()

    save_pred_results = para_dict['save_folder']+'/pred_results.pkl'
    print('The prediction results have been generated and saved under: ', save_pred_results)
    save_file(save_pred_results, pred_results)



    '''
    Draw the ROC curves .... 
    '''
    # fpr, tpr, thresholds = roc_curve(y_true_list, y_score_list)
    # AUC_score = roc_auc_score(y_true_list, y_score_list)
    # save_ROC_results = para_dict['save_folder']+'/'+para_dict['save_folder'].split('/')[len(para_dict['save_folder'].split('/'))-2]+'__ROC_results.pkl'
    # save_file(save_ROC_results, (fpr, tpr, thresholds, AUC_score))
    # print('XXXXXXX The AUC score is: ', round(AUC_score, 2))
    # fig, ax = plt.subplots(figsize=(20, 10))
    # plt.plot(fpr, tpr, 'r', linewidth=2, label='ROC curve [Area Under Curve (AUC = {:.4f})]'.format(AUC_score))
    # plt.legend(loc='lower right', fontsize=20)
    # plt.plot([0, 1], [0, 1], color='navy', linestyle='--')
    # plt.xlim([-0.02, 1.02])
    # plt.ylim([-0.02, 1.02])
    # plt.xlabel('Specificity', fontsize=20)
    # plt.ylabel('Sensitivity', fontsize=20)
    # # plt.xlabel('False Positive Rate', fontsize=20)
    # # plt.ylabel('True Positive Rate', fontsize=20)
    # n_split = 11
    # xTick_list = []
    # for n in np.linspace(0, 1, n_split):
    #     xTick_list.append(str(int(n*100))+'%')
    # # reversing the list
    # new_xTick_list = []
    # for i in xTick_list:
    #     new_xTick_list.insert(0, i)
    # plt.xticks(np.linspace(0, 1, n_split), new_xTick_list, fontsize=15)
    # yTick_list = []
    # for n in np.linspace(0, 1, n_split):
    #     yTick_list.append(str(int(n*100))+'%')
    # plt.yticks(np.linspace(0, 1, n_split), yTick_list, fontsize=15)
    # plt.grid(color='y', linewidth=0.5)
    # # plt.title('Mean ROC curve for TB Index Score', fontsize=35)
    # plt.show()
    # plt.savefig(saveFig_READ_ROC_title)
    # print('The ROC figure has been saved at: ', saveFig_READ_ROC_title)
    # plt.close('all')


    for mapc in pred_results:
        print("Total sites for <{}>: +{} -{}\n".format( mapc, (len(pred_results[mapc]['+']) if '+' in pred_results[mapc] else 0), (len(pred_results[mapc]['-']) if '-' in pred_results[mapc] else 0) ) );
        for mapstn in pred_results[mapc]:
            for ref_pos in pred_results[mapc][mapstn]:
                
                t_pred = pred_results[mapc][mapstn][ ref_pos ]
                # if ref_pos == 2161:
                #     print('t_pred: ', t_pred)
                
                print("Pred: {}:{}:{} {} {} {} {}".format( mapc, mapstn, ref_pos, t_pred[0][0], t_pred[0][1], t_pred[1][0], t_pred[1][1] ) )
    print()
    read_pred_writer.close()
    sys.stdout.flush()

    coverage_thr = 1;
    #             pos      neg
    #            tp fn    fp tn            
    site_perf = [[0, 0], [0, 0]]
    site_perf_detail = [[], []]
    read_perf = [[0, 0], [0, 0]]
    perc_thr = 0.35;

    for mapc in pred_results:
        for mapstrnd in pred_results[mapc]:
            for mappos in pred_results[mapc][mapstrnd]:
                t_pred_info = pred_results[ mapc ][ mapstrnd ][ mappos ]
                for _pnind in range(2):
                    # print('t_pred_info[_pnind][0] + t_pred_info[_pnind][1]: ', t_pred_info[_pnind][0] + t_pred_info[_pnind][1])
                    if t_pred_info[_pnind][0] + t_pred_info[_pnind][1] >= coverage_thr:
                        read_perf[_pnind][0] += t_pred_info[_pnind][0]
                        read_perf[_pnind][1] += t_pred_info[_pnind][1]
                for _pnind in range(2):
                    if t_pred_info[_pnind][0] + t_pred_info[_pnind][1] >= coverage_thr:
                        site_perf_detail[ _pnind ].append( t_pred_info[_pnind][0]/float(t_pred_info[_pnind][0] + t_pred_info[_pnind][1] ) )
                        if t_pred_info[_pnind][0]/float(t_pred_info[_pnind][0] + t_pred_info[_pnind][1]) < perc_thr:
                            site_perf[ _pnind ][1] += 1;
                        else:
                            site_perf[ _pnind ][0] += 1;
    
    print("site_perf: {}".format(site_perf))
    calculate_perf(site_perf)
    print("read_perf: {}".format(read_perf))
    calculate_perf( read_perf )


