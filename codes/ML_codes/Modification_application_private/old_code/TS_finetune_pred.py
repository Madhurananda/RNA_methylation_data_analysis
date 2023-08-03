
import os,sys
import random
import time, datetime
import numpy as np
import copy
import math
import traceback


import torch
import torch.multiprocessing


sys.path.append( os.path.dirname(__file__) )

#import TS_model
import TS_biLSTM
#import TS_dataset_reader

import TS_global
import TS_ftdataset_reader
import read_tombo


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

        a_pos_list = []
        r_pos_list = []
        for posi in range(map_end - map_sta):
            if map_stn=='+':
                #if ref_seq_dict[ map_chr ][ posi + map_sta ] in ['T', 't']: #['A', 'a']:
                if ref_seq_dict[ map_chr ][ posi + map_sta ] in ['A', 'a']:
                    a_pos_list.append( posi )
                    r_pos_list.append( posi + map_sta )
                    #print("test: {} {} {}/{} {} {} ".format( map_end, map_sta, map_stn,  posi, map_chr, posi + map_sta ))
            else:
                #if ref_seq_dict[ map_chr ][ map_end - posi - 1] in ['A', 'a']: # ['T', 't']:
                if ref_seq_dict[ map_chr ][ map_end - posi - 1] in ['T', 't']:
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
    print("All chrs={}".format( sorted(list(seq_dict.keys())) ) );
    print("Use chrs={}".format( use_chr ) );
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
    used_cpu_devices = 3

    start_time = time.time()
    pred_results = {}
    mod_signal_model.eval();
    loss_total = [0, 0]
    pred_size = [0, 0]

    read_pred_writer = open(para_dict['save_file'], 'w')

    process_i = 0;
    with torch.no_grad(): 
        for dataf_ind in range(len(datafolders)):
            dataf = datafolders[ dataf_ind ]
            fs_dataset_0 = TS_ftdataset_reader.TS_Dataset(dataf, index_file)
            fs_dataloader_0 = torch.utils.data.DataLoader(fs_dataset_0, batch_size=1, num_workers=used_cpu_devices, shuffle=True)
            print(type(fs_dataset_0), len(fs_dataset_0), fs_dataset_0.total_file_seq_count, fs_dataset_0.total_seq_count )
            sys.stdout.flush()

            testi = 0;
            for batch_fidx_0, batch_data_0 in enumerate( fs_dataloader_0 ):
                feat_ts_dict, pad_inf_dict, map_info_dict, fn_list = generate_pred_feat ( batch_data_0, seq_dict, seqr, para_dict)
                #if process_i<1 and testi<10:
                #    print("{}/{} {} {} {}".format(testi, process_i ,len(fn_list), len(feat_ts_dict), len(map_info_dict) ))
                #    sys.stdout.flush()
                #sys.exit(0)
                if len( fn_list ) ==0: continue;
                
                testi += 1; 
                #if testi >= fs_dataset_0.total_file_seq_count*0.1: break;

                for _fn_ in fn_list:
                    process_i += 1;
                    map_info = map_info_dict[_fn_]
                    if map_info[0] not in pred_results:
                        pred_results[ map_info[0] ] = {};
                        #print('add {} / {}'.format( map_info[0], pred_results.keys() ) )
                    if map_info[1] not in pred_results[ map_info[0] ]:
                        pred_results[ map_info[0] ][ map_info[1] ] = {};
                        #print('add {} / {}'.format( map_info[1], pred_results[ map_info[0] ].keys()) )
                    #        1  2  3  4             1        2        3            4
                    read_pred_writer.write(">{} {} {} {}\n".format( dataf_ind, _fn_, map_info[0], map_info[1] ))

                    #print ( type(map_info[4]), type(feat_ts_dict[_fn_]), map_info[:2] )
                    #print(torch.stack(feat_ts_dict[_fn_], 0).size())
                    #print (len(feat_ts_dict[_fn_]))
                    #print ( len(map_info[4]), map_info[:2] )
                    #read_pred_writer.close()
                    #sys.exit(0)
                    #<class 'list'> <class 'list'> ('chrV', '-')
                    #torch.Size([1976, 51, 100])
                    #1976
                    #1976 ('chrV', '-')


                    ###
                    for start_ind in range(0, len(feat_ts_dict[_fn_]), input_batch_size): 
                        end_ind = len(feat_ts_dict[_fn_]) if start_ind+input_batch_size*1.5>=len(feat_ts_dict[_fn_]) else start_ind+input_batch_size
                        t_feat = torch.stack(feat_ts_dict[_fn_][start_ind:end_ind], 0)
                        pred_mask = mod_signal_model.forward( t_feat.to(device), device );
                        if dataf_ind==0:
                            class_lab = torch.zeros(t_feat.size(0), dtype=torch.long ).to(device)
                        else: 
                            class_lab = torch.ones( t_feat.size(0), dtype=torch.long ).to(device)

                        #print(_fn_, start_ind, len(feat_ts_dict[_fn_]), pred_mask.size(), class_lab.size() )
                        loss = loss_func(pred_mask, class_lab)
                        loss_total[ dataf_ind ] += loss*t_feat.size(0)
                        pred_size[ dataf_ind ] += t_feat.size(0)

                        pred_id = torch.max(pred_mask, 1)[1]
                        #for ref_pos_ind in range(len(map_info[4])):
                        for ref_pos_ind in range(start_ind, end_ind):    
                            ref_pos = map_info[4][ ref_pos_ind ]
                            read_pred_writer.write("{} {}\n".format( ref_pos, pred_id[ref_pos_ind-start_ind].item() ))
                            if ref_pos not in pred_results[ map_info[0] ][ map_info[1] ]:
                                pred_results[ map_info[0] ][ map_info[1] ][ ref_pos ] = [[0, 0], [0, 0]]
                            pred_results[ map_info[0] ][ map_info[1] ][ ref_pos ][ dataf_ind ] [ pred_id[ref_pos_ind-start_ind].item() ] += 1
                        if end_ind==len(feat_ts_dict[_fn_]):
                            break;
                    ###
                    read_pred_writer.write(">\n")

                    if process_i%5000==0:
                        print("{} Processing {} for {:.1f} {:.5f}/{:,} {}".format(dataf_ind, process_i, time.time()-start_time, loss_total[ dataf_ind ]/(pred_size[ dataf_ind ] if pred_size[ dataf_ind ]>0 else 1), pred_size[ dataf_ind ], t_feat.shape))
                        start_time = time.time()
                        sys.stdout.flush()
                    #break;
                #break;
            print("{} Processing {} for {:.1f} {:.5f}/{:,} ".format(dataf_ind, process_i, time.time()-start_time, loss_total[ dataf_ind ]/(pred_size[ dataf_ind ] if pred_size[ dataf_ind ]>0 else 1), pred_size[ dataf_ind ] ))
            start_time = time.time()
            sys.stdout.flush()
    for mapc in pred_results:
        print("Total sites for <{}>: +{} -{}\n".format( mapc, (len(pred_results[mapc]['+']) if '+' in pred_results[mapc] else 0), (len(pred_results[mapc]['-']) if '-' in pred_results[mapc] else 0) ) );
        for mapstn in pred_results[mapc]:
            for ref_pos in pred_results[mapc][mapstn]:
                t_pred = pred_results[mapc][mapstn][ ref_pos ]
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
                t_pred_info = pred_results[mapc][mapstrnd][ mappos ]
                for _pnind in range(2):
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


#if __name__=='__main__':
if __name__=='__main__' and '__file__' in globals():
    para_dict = {}
    para_dict['save_folder'] = './testPred';
    if not os.path.isdir( para_dict['save_folder'] ):
        os.system( 'mkdir -p {}'.format( para_dict['save_folder'] ) )
    #para_dict['save_file_pref'] = 'BilstmMean.'+'layers'+sys.argv[1]+'.hs'+sys.argv[2]+'.'+sys.argv[4]+'.lr'+sys.argv[5]+'.b'+sys.argv[6]+'.'
    basefolder = '/home2/qliu10/Transignaler_bp';

    para_dict['size_layers'] =int(sys.argv[1])
    para_dict['size_hidden'] =int(sys.argv[2])

    para_dict['cpu_devices'] = 15; #
    para_dict['cuda_devices'] = 1;
    para_dict['cuda_id'] = sys.argv[3];

    saved_model = sys.argv[5]
    print(datetime.datetime.now())
    try:
        datafs = [basefolder+'/tfrecord_ft.5M/Curlcake_m6a_cc6m_2244_T7_ecorv/', basefolder+'/tfrecord_ft.5M/Curlcake_non_cc6m_2244_T7_ecorv/']
        pred( para_dict = para_dict, datafolders=datafs, saved_model=saved_model, ref_seq_file='EpiNano/Reference_sequences/cc.fasta',index_file="tfrecord.index", input_batch_size=int(sys.argv[4]), use_chr=['cc6m_2244_T7_ecorv'] )
    except Exception as e:
        print(traceback.format_exc())
        print(e)
        print(datetime.datetime.now())
    print(datetime.datetime.now())


