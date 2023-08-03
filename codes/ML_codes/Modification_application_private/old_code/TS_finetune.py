
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

##################
###
def generate_train_feat(signal_full_info2, p_batch_size, used_cpu_devices, ref_seq_dict, seqr, para_dict):
    feat_ts = []
    pad_inf = []

    total_signal_event = 0
    signal_full_length = copy.deepcopy(signal_full_info2[0])
    for _i_seq in range(len(signal_full_length)):
        if signal_full_length[_i_seq].size(1)<para_dict['length_thr']: continue;
        map_chr = signal_full_info2[3][_i_seq][0][0]
        if map_chr not in ref_seq_dict:
            continue;
        map_stn = signal_full_info2[3][_i_seq][1][0]
        map_sta = signal_full_info2[3][_i_seq][2]
        map_end = signal_full_info2[3][_i_seq][3]

        ####
        ### pos highly rely on the motif; 
        ### need revise
        ###
        a_pos_list = []
        r_pos_list = []
        for posi in range(map_end - map_sta):
            if map_stn=='+':
                #if ref_seq_dict[ map_chr ][ posi + map_sta ] in ['T', 't']: # ['A', 'a']:
                if ref_seq_dict[ map_chr ][ posi + map_sta ] in ['A', 'a']:
                    a_pos_list.append( posi )
                    r_pos_list.append( posi + map_sta )
            else:
                #if ref_seq_dict[ map_chr ][ map_end - posi - 1] in ['A', 'a']: # ['T', 't']:
                if ref_seq_dict[ map_chr ][ map_end - posi - 1] in ['T', 't']:
                    a_pos_list.append( posi )
                    r_pos_list.append( map_end - posi - 1 )

        total_signal_event += signal_full_length[_i_seq].size(1)
        #if signal_full_length[_i_seq].size(1)<TS_global.min_seq_len: continue;
        #if signal_full_length[_i_seq].size(1)<para_dict['length_thr']: continue;

        #random_center_pos = np.random.randint(1, TS_global.input_len//2);
        #for _this_cp in range(random_center_pos, signal_full_length[_i_seq].size(1), random_center_pos):
        for _this_cp in a_pos_list:
            if _this_cp<5 or _this_cp>signal_full_length[_i_seq].size(1)-5: continue;

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

    pad_inf = torch.tensor(pad_inf)
    feat_ts = torch.stack(feat_ts)

    train_dataset = torch.utils.data.TensorDataset(feat_ts, pad_inf)
    train_loader = torch.utils.data.DataLoader(train_dataset, batch_size=p_batch_size, collate_fn=train_feat_Batch_wrapper, pin_memory=True, shuffle=True, num_workers=used_cpu_devices)

    return train_loader;


def get_adapt_learning_rate(init_lr, warmup_seq, num_train_seq, max_estimat_seq):
    decay_rate = 1;
    if num_train_seq <= warmup_seq:
        decay_rate = num_train_seq/warmup_seq
    else:
        #if num_train_seq>max_estimat_seq-max_estimat_seq*0.01:
        if num_train_seq-warmup_seq>(max_estimat_seq-warmup_seq)*0.99:
            #decay_rate = np.min([1- max_estimat_seq*0.01/(max_estimat_seq-warmup_seq),
            #                 np.power( (max_estimat_seq*0.01/(warmup_seq)), -0.5 )]);
            decay_rate = math.cos( math.pi/2*(max_estimat_seq-warmup_seq)*0.99/(max_estimat_seq-warmup_seq) )
        else:
            #decay_rate = np.min([1- ( num_train_seq-warmup_seq )/(max_estimat_seq-warmup_seq),
            #                np.power( ((num_train_seq-warmup_seq )/(warmup_seq)), -0.5 )]);
            decay_rate = math.cos( ( math.pi/2*(num_train_seq-warmup_seq) )/(max_estimat_seq-warmup_seq) )
    return decay_rate*init_lr

def save_mod_model(mod_signal_model, _optmizer, device, para_dict, epoch, num_train_seq, loss_func, num_train_1epoch, max_estimat_seq):
    output_path = para_dict['save_folder']+'/'+para_dict['save_file_pref'] + ".GPU.ep{}.{}.pt".format(epoch, num_train_seq)
    torch.save({
                "epoch": epoch,
                "num_train_seq": num_train_seq,
                "model_state_dict": mod_signal_model.state_dict(),
                "optimizer_state_dict": _optmizer.state_dict(),
                "loss_func": loss_func,
                "num_train_1epoch": num_train_1epoch,
                "max_estimat_seq": max_estimat_seq
                }, 
                output_path)
    if 'cpumodel' in para_dict:
        output_path = para_dict['save_folder']+'/'+para_dict['save_file_pref'] + ".CPU.ep{}.{}.pt".format(epoch, num_train_seq)
        torch.save({
                    "epoch": epoch,
                    "num_train_seq": num_train_seq,
                    "model_state_dict": mod_signal_model.cpu().state_dict(), 
                    "optimizer_state_dict": _optmizer.state_dict(),
                    "loss_func": loss_func,
                    "num_train_1epoch": num_train_1epoch,
                    "max_estimat_seq": max_estimat_seq
                    }, 
                    output_path)
        mod_signal_model.to(device)
        #_optmizer.to(device)
    print("EP:{}/{} Model Saved on: {}".format( epoch, num_train_seq, output_path));


def finetuneTrainer(para_dict, datafolders, is_ft, ref_seq_file, not_use_chr=[], index_file='tfrecord.index', input_batch_size=32, cuda_devices: int=1, cpu_devices: int=5, adam_learning_rate: float =1e-4, ft_learning_rate=1e-6, adam_betas=(0.9,0.999), adam_weight_decay: float=0.01, random_seed = 1,  train_epoch=50, warmup_seq=10000, saved_model=None, downsampling_rate=None):
    torch.manual_seed(random_seed)
    np.random.seed(random_seed)
    random.seed(random_seed)

    torch.multiprocessing.set_sharing_strategy('file_system')

    if not os.path.isdir(para_dict['save_folder']):
        os.system('mkdir -p {}'.format(para_dict['save_folder']))

    #seqstep = para_dict['save_frequency']
    #seqstep_small = para_dict['save_frequency_s']
    seqstep = para_dict['save_frequency']
    logstep = para_dict['log_frequency']
    if 'saved_model' in para_dict:
        saved_model = para_dict['saved_model']
    if 'warmup_seq' in para_dict:    
        if para_dict['warmup_seq']>100: warmup_seq    = para_dict['warmup_seq']
        else: print('warmup_seq is too small: {}. not used.'.format(para_dict['warmup_seq']))
    if 'train_epoch' in para_dict:   
        if para_dict['train_epoch']>0: train_epoch   = para_dict['train_epoch']
        else: print('train_epoch is too small: {}. not used.'.format(para_dict['train_epoch']))
    if 'random_seed' in para_dict:   
        random_seed   = para_dict['random_seed']
    if 'cpu_devices' in para_dict:   
        if para_dict['cpu_devices']>0: cpu_devices   = para_dict['cpu_devices']
        else: print('cpu_devicesis too small: {}. not used.'.format(para_dict['cpu_devices']))
    if 'cuda_devices' in para_dict:  
        if para_dict['cuda_devices']>0: cuda_devices  = para_dict['cuda_devices']
        else:  print('cuda_devices too small: {}. not used.'.format(para_dict['cuda_devices']))
    if 'input_batch_size' in para_dict: 
        if para_dict['input_batch_size']>0: input_batch_size = para_dict['input_batch_size']
        else: print('input_batch_sizeis too small: {}. not used.'.format(para_dict['input_batch_size']))
    if 'index_file' in para_dict:    
        index_file    = para_dict['index_file']
    if 'ft_learning_rate' in para_dict:
        ft_learning_rate = para_dict['ft_learning_rate']

    if 'adam_learning_rate' in para_dict:
        adam_learning_rate = para_dict['adam_learning_rate']
    if 'downsampling_rate' in para_dict:
        downsampling_rate = para_dict['downsampling_rate']
    
    input_batch_size2 = input_batch_size

    #print("Key learn info input_batch_size={} adam_learning_rate={} warmup_seq={:,} seqstep={:,} seqstep_small={:,}".format(input_batch_size, adam_learning_rate, warmup_seq, seqstep, seqstep_small))
    print("Key learn info input_batch_size={} adam_learning_rate={} warmup_seq={:,} seqstep={:,} logstep={:,} cudaID={} size_layers={} size_hidden={} is_ft={}/{} {}".format(input_batch_size, adam_learning_rate, warmup_seq, seqstep, logstep, ("cuda:"+para_dict['cuda_id'] if 'cuda_id' in para_dict else "cuda:0"), (3 if 'size_layers' not in para_dict or para_dict['size_layers']<1 else para_dict['size_layers']), (1024 if 'size_hidden' not in para_dict or para_dict['size_hidden']<10 else para_dict['size_hidden']), is_ft, (ft_learning_rate), ("" if not is_ft else saved_model) ))

    if is_ft and saved_model==None:
        print("Error!!! No model for fine-tune")
        return;

    #mod_signal_model = TS_model.SignalBertPred( );
    mod_signal_model = TS_biLSTM.SignalPred( ispretrain = False, pretrainmean=False,
                                                 size_layers = (3 if 'size_layers' not in para_dict or para_dict['size_layers']<1 else para_dict['size_layers']),
                                                 size_hidden = (1024 if 'size_hidden' not in para_dict or para_dict['size_hidden']<10 else para_dict['size_hidden'] ) );
    #device = torch.device("cuda:0" if torch.cuda.is_available() else "CPU")
    device = torch.device(("cuda:"+para_dict['cuda_id'] if 'cuda_id' in para_dict else "cuda:0") if torch.cuda.is_available() else "CPU")
    mod_signal_model = mod_signal_model.to(device);
    '''para_sum = 0;
    pretrain_param = []
    ft_param = []
    pretrain_param_names = []
    ft_param_names = []
    for l_name, l_para in mod_signal_model.named_parameters():
        print(l_name, l_para.data.shape, l_para.requires_grad,  type(l_name) )
        if l_name.split('.')[0] not in ['mask_linear', 'mask_linear1', 'mask_linear2', 'mask_linear1_2']:
            para_sum += l_para.nelement(); #l_para.data.size(0);
        if 'pred_' in l_name:
            ft_param.append(l_para)
            ft_param_names.append(l_name)
        else:
            pretrain_param.append(l_para)
            pretrain_param_names.append(l_name)
    print("Total parameters: {}".format( para_sum ))   
    print(ft_param_names, len(ft_param)) ; #, ft_param)
    print('\n')
    print(pretrain_param_names, len(pretrain_param)); #, pretrain_param)
    print('\n')
    '''
    cpuCount = os.cpu_count()
    used_cpu_devices = cpu_devices
    if cpuCount <= used_cpu_devices : used_cpu_devices = cpuCount
    print("Number of CPUs in the system: {}. Used CPUs: {}".format(cpuCount, used_cpu_devices ) )

    '''if torch.cuda.device_count() > 1 and cuda_devices>0:
        if cuda_devices>=torch.cuda.device_count(): used_cuda_devices = torch.cuda.device_count()
        else: used_cuda_devices = cuda_devices
        print("Number of GPUs in the system: {}. {} GPUs for training".format( torch.cuda.device_count(), used_cuda_devices ))
        mod_signal_model = torch.nn.DataParallel( mod_signal_model, device_ids=list(range( used_cuda_devices )) );
    '''
    max_estimat_seq = 0;
    num_train_seq = 0;
    num_train_i = 0;
    save_model_f = 0;
    first_n_warmup = 3;
    train_epoch_i = 0
    num_train_1epoch = 0;

    seq_dict = TS_ftdataset_reader.read_seq_id( ref_seq_file )
    print("All chrs={}".format( sorted(list(seq_dict.keys())) ) );
    print("Nus chrs={}".format( not_use_chr ) );
    for nous_c in not_use_chr:
        if nous_c not in seq_dict:
            print("Warning!!! Specified no-used Chr <{}> is not in ref files={}".format( nous_c, ref_seq_file ))
        else:
            del seq_dict[ nous_c ]
    seqr = read_tombo.SeqReverse()

    if not saved_model==None:
        t_checkpoint = torch.load( saved_model ); #, map_location=device );
        mod_signal_model.load_state_dict(t_checkpoint['model_state_dict'])
        '''_optmizer.load_state_dict(t_checkpoint['optimizer_state_dict'])
        num_train_seq = t_checkpoint['num_train_seq']
        save_epoch_i = t_checkpoint['epoch']
        train_epoch_i = t_checkpoint['epoch']-1
        loss_func = t_checkpoint['loss_func']
        num_train_1epoch = t_checkpoint['num_train_1epoch']
        max_estimat_seq = t_checkpoint['max_estimat_seq']
        total_max_estimat_seq = max_estimat_seq

        load_num_train_seq = 0;
        for _prev_epoch_i in range(0, train_epoch_i):
            load_num_train_seq += num_train_1epoch
        '''
    para_sum = 0;
    pretrain_param = []
    ft_param = []
    pretrain_param_names = []
    ft_param_names = []
    for l_name, l_para in mod_signal_model.named_parameters():
        print(l_name, l_para.data.shape, l_para.requires_grad,  type(l_name) )
        if l_name.split('.')[0] not in ['mask_linear', 'mask_linear1', 'mask_linear2', 'mask_linear1_2']:
            para_sum += l_para.nelement(); #l_para.data.size(0);
        if 'pred_' in l_name:
            ft_param.append(l_para)
            ft_param_names.append(l_name)
        else:
            pretrain_param.append(l_para)
            pretrain_param_names.append(l_name)
    print("Total parameters: {}".format( para_sum ))
    print(ft_param_names, len(ft_param)) ; #, ft_param)
    print('\n')
    print(pretrain_param_names, len(pretrain_param)); #, pretrain_param)
    print('\n')

    if is_ft:
        _optmizer = torch.optim.Adam( [{'params':pretrain_param}, {'params':ft_param}], lr = adam_learning_rate, betas=adam_betas, weight_decay=adam_weight_decay )
        _optmizer.param_groups[0]['lr'] = ft_learning_rate
        _optmizer.param_groups[1]['lr'] = adam_learning_rate
        if True:
            _optmizer = torch.optim.SGD(  [{'params':pretrain_param}, {'params':ft_param}], lr = adam_learning_rate, momentum=0.9);
            _optmizer.param_groups[0]['lr'] = ft_learning_rate
            _optmizer.param_groups[1]['lr'] = adam_learning_rate
    else: 
        #_optmizer = torch.optim.Adam( mod_signal_model.parameters(), lr = adam_learning_rate, betas=adam_betas, weight_decay=adam_weight_decay )
        _optmizer = torch.optim.SGD( mod_signal_model.parameters(), lr = adam_learning_rate, momentum=0.9);
    loss_func = torch.nn.CrossEntropyLoss()

    used_cpu_devices = 0

    #try no fine-tune
    #for param_t in _optmizer.param_groups[0]['params']:
    #    param_t.requires_grad = False;

    sum_los = 0;
    acc_val = 0;
    previous_num_train_i = num_train_i
    start_time = time.time()
    start_adam_learning_rate = adam_learning_rate;
    #for train_epoch_i in range(train_epoch):
    fs_dataset_1 = None;
    while train_epoch_i < train_epoch:
        train_epoch_i += 1; # epoch_i starts from 1 to train_epoch
        #adam_learning_rate = start_adam_learning_rate/(2**(train_epoch_i-1))
        adam_learning_rate = start_adam_learning_rate/((np.sqrt(2)**(train_epoch_i-1)) if train_epoch<10 else np.sqrt(train_epoch_i) )

        fs_dataset_0 = TS_ftdataset_reader.TS_Dataset(datafolders[0], index_file, m_seed=train_epoch_i*random_seed, use_percentage=downsampling_rate)
        if train_epoch_i==1:
           #if fs_dataset_0.total_file_seq_count>=len(fs_dataset_0)-10 and fs_dataset_0.total_file_seq_count>1 :
           #    max_estimat_seq = fs_dataset_0.total_seq_count
           #else:
           #    max_estimat_seq = len(fs_dataset_0)*(800 if TS_global.tf_record_max_event ==1000000 else 400)
           max_estimat_seq = fs_dataset_0.total_seq_count
           total_max_estimat_seq = max_estimat_seq
           print(type(fs_dataset_0), len(fs_dataset_0), fs_dataset_0.total_file_seq_count, fs_dataset_0.total_seq_count )
           sys.stdout.flush()
        fs_dataloader_0 = torch.utils.data.DataLoader(fs_dataset_0, batch_size=1, num_workers=used_cpu_devices, shuffle=True)

        #fs_dataset_1 = TS_ftdataset_reader.TS_Dataset(datafolders[1], index_file, m_seed=train_epoch_i*random_seed)
        #fs_dataloader_1 = None;
        #fs_dataloader_1 = torch.utils.data.DataLoader(fs_dataset_1, batch_size=1, num_workers=used_cpu_devices, shuffle=True)

        print("Epoch: {}/{} new_learning_rate={:.1e}".format( train_epoch_i, train_epoch, adam_learning_rate))
        for batch_fidx_0, batch_data_0 in enumerate( fs_dataloader_0 ):
            t_seq_num = len(batch_data_0[0])
            #if (not saved_model==None) and load_num_train_seq+t_seq_num<=num_train_seq:
            #    load_num_train_seq +=t_seq_num
            #    continue;
            save_model_f = 0;
            if (num_train_seq+t_seq_num)//logstep> (num_train_seq//logstep):
                save_model_f = 1;
            if ((num_train_seq+t_seq_num)//(seqstep)) > (num_train_seq//(seqstep)):
                save_model_f = 2
            num_train_seq += t_seq_num
            if train_epoch_i==1:
                num_train_1epoch += t_seq_num
            #new_lr = get_adapt_learning_rate(adam_learning_rate, warmup_seq, num_train_seq-total_max_estimat_seq*(train_epoch_i-1), total_max_estimat_seq);
            new_lr = adam_learning_rate; #get_adapt_learning_rate(adam_learning_rate, warmup_seq, num_train_seq-total_max_estimat_seq*(train_epoch_i-1), total_max_estimat_seq);
            if is_ft:
               if num_train_i<1:
                   print( _optmizer )
                   print()
                   print(_optmizer.param_groups[0] )
                   print()
                   print( _optmizer.param_groups[1] )
                   sys.stdout.flush()
               #if not is_ft: _optmizer.param_groups[0]['lr'] = new_lr
               _optmizer.param_groups[1]['lr'] = new_lr
               if train_epoch_i>train_epoch//2:
                   for param_t in _optmizer.param_groups[0]['params']:
                        pass; #param_t.requires_grad = False;
            else:
                for param_group in _optmizer.param_groups:
                    param_group['lr'] = new_lr

            train_loader_0 = generate_train_feat( batch_data_0, input_batch_size, used_cpu_devices, seq_dict, seqr, para_dict)
            ######################################
            #if fs_dataloader_1==None:
            #    fs_dataloader_1 = torch.utils.data.DataLoader(fs_dataset_1, batch_size=1, num_workers=used_cpu_devices, shuffle=True)
            #batch_data_1 = next(iter(fs_dataloader_1))
            #if not batch_data_1:
            #    fs_dataloader_1 = torch.utils.data.DataLoader(fs_dataset_1, batch_size=1, num_workers=used_cpu_devices, shuffle=True)
            #    batch_data_1 = next(iter(fs_dataloader_1))
            #    train_loader_1 = generate_train_feat( batch_data_1, input_batch_size2, used_cpu_devices, seq_dict, seqr))
            ######################################

            for batch_ndx_0, train_batch_0 in enumerate(train_loader_0):
                if train_batch_0.feat.size(0) < input_batch_size/2: continue;
                num_train_i += 1
        
                if fs_dataset_1==None:
                    fs_dataset_1 = TS_ftdataset_reader.TS_Dataset(datafolders[1], index_file, m_seed=train_epoch_i*random_seed, use_percentage=downsampling_rate)
                    fs_dataloader_1 = torch.utils.data.DataLoader(fs_dataset_1, batch_size=1, num_workers=used_cpu_devices, shuffle=True)

                    batch_data_1 = next(iter(fs_dataloader_1))
                    train_loader_1 = generate_train_feat( batch_data_1, input_batch_size2, used_cpu_devices, seq_dict, seqr, para_dict)
                train_batch_1 =  next(iter(train_loader_1))
                while (train_batch_1) and train_batch_1.feat.size(0) < input_batch_size2/2:
                    train_batch_1 =  next(iter(train_loader_1))
                while (not train_batch_1):
                    batch_data_1 = next(iter(fs_dataloader_1))
                    if not batch_data_1:
                        fs_dataset_1 = TS_ftdataset_reader.TS_Dataset(datafolders[1], index_file, m_seed=train_epoch_i*random_seed, use_percentage=downsampling_rate)
                        fs_dataloader_1 = torch.utils.data.DataLoader(fs_dataset_1, batch_size=1, num_workers=used_cpu_devices, shuffle=True)

                        batch_data_1 = next(iter(fs_dataloader_1))
                    train_loader_1 = generate_train_feat( batch_data_1, input_batch_size2, used_cpu_devices, seq_dict, seqr, para_dict)
                    train_batch_1=  next(iter(train_loader_1))

                #####################
                _optmizer.zero_grad();

                class_lab_0 = torch.zeros( train_batch_0.feat.size(0), dtype=torch.long ) 
                class_lab_1 = torch.ones( train_batch_1.feat.size(0), dtype=torch.long )
                class_lab = torch.cat( (class_lab_0, class_lab_1), 0 ).to(device)
                tb_feat = torch.cat( (train_batch_0.feat, train_batch_1.feat), 0).to(device)
                pred_mask = mod_signal_model.forward(tb_feat, device)
                #print( pred_mask.shape, class_lab.shape)
                loss = loss_func(pred_mask, class_lab) 
    
                #_optmizer.zero_grad(); 
                loss.backward();
                _optmizer.step();

                if num_train_i in [5000, 20000, 50000, 80000, 200000, 500000, 800000]:
                    print( pred_mask )

                sum_los += loss.item()
                acc_val += (torch.max(pred_mask, 1)[1]==class_lab).sum().item()/float( pred_mask.size(0))
                if num_train_i<100 and num_train_i-previous_num_train_i>3:
                    #print("This loss: {:.3f} for ep{}.{:,}. Consuming time: {:.1f}. Learning_rate={:.5e} avgloss={:.3f} traini={:,} acc={} {}/{} {}/{} {}".format(loss.item(), train_epoch_i, num_train_seq, time.time()-start_time, new_lr, sum_los/(num_train_i-previous_num_train_i), num_train_i , acc_val/(num_train_i-previous_num_train_i) , pred_mask.shape, pred_mask, class_lab.shape, class_lab, loss.item()))
                    print("This loss: {:.3f} for ep{}.{:,}. Consuming time: {:.1f}. Learning_rate={:.5e} avgloss={:.3f} traini={:,} acc={:.3f} {} {} {}".format(loss.item(), train_epoch_i, num_train_seq, time.time()-start_time, new_lr, sum_los/(num_train_i-previous_num_train_i), num_train_i , acc_val/(num_train_i-previous_num_train_i) , pred_mask.shape,  class_lab.shape, loss.item()))

            # end for batch_ndx, train_batch in enumerate(train_loader):
            if (save_model_f>0 or (num_train_seq<seqstep//10) ) and num_train_i-previous_num_train_i>10:
                #print("This loss: {} for ep{}.{}. Consuming time: {:.1f}".format(loss.item(), train_epoch_i, num_train_seq, time.time()-start_time))
                #start_time = time.time()
                #save_mod_model(mask_signal_model, _optmizer, device, para_dict, train_epoch_i, num_train_seq, loss_func)
                print("This loss:: {:.3f} for ep{}.{:,}. Consuming time: {:.1f}. Learning_rate={:.5e} avgloss={:.3f} traini={:,} acc={:.3f}".format(loss.item(), train_epoch_i, num_train_seq, time.time()-start_time, new_lr, sum_los/(num_train_i-previous_num_train_i), num_train_i, acc_val/(num_train_i-previous_num_train_i) ))
                start_time = time.time()
                sum_los = 0;
                acc_val = 0
                previous_num_train_i = num_train_i
                if save_model_f>1:
                    save_mod_model(mod_signal_model, _optmizer, device, para_dict, train_epoch_i, num_train_seq, loss_func, num_train_1epoch, max_estimat_seq)
                    #save_mod_model(mod_signal_model, _optmizer, device, para_dict,       epoch,   num_train_seq, loss_func, num_train_1epoch, max_estimat_seq)
                sys.stdout.flush()
        if not ( (save_model_f>0 or (num_train_seq<seqstep//10) ) and num_train_i-previous_num_train_i>10 ):
            if not (previous_num_train_i==num_train_i):
                print("This loss:: {:.3f} for ep{}.{:,}. Consuming time: {:.1f}. Learning_rate={:.5e} avgloss={:.3f} traini={:,} acc={:.3f}".format(loss.item(), train_epoch_i, num_train_seq, time.time()-start_time, new_lr, sum_los/(num_train_i-previous_num_train_i), num_train_i, acc_val/(num_train_i-previous_num_train_i) ))
                start_time = time.time()
                sum_los = 0;
                acc_val = 0
                previous_num_train_i = num_train_i
        save_mod_model(mod_signal_model, _optmizer, device, para_dict, train_epoch_i, num_train_seq, loss_func, num_train_1epoch, max_estimat_seq)
        #end for batch_fidx, batch_data in enumerate( fs_dataloader ):
    #save_mod_model(mask_signal_model, _optmizer, device, para_dict, train_epoch_i, num_train_seq, loss_func)
    #save_mod_model(mod_signal_model, _optmizer, device, para_dict, train_epoch_i, num_train_seq, loss_func, num_train_1epoch, max_estimat_seq)

#if __name__=='__main__':
if __name__=='__main__' and '__file__' in globals():
    para_dict = {}
    para_dict['save_folder'] = './testModel';
    para_dict['save_file_pref'] = 'BilstmMean.'+'layers'+sys.argv[1]+'.hs'+sys.argv[2]+'.'+sys.argv[4]+'.lr'+sys.argv[5]+'.b'+sys.argv[6]+(".p"+str(int(float(sys.argv[7])*1000+0.5)) if len(sys.argv)>7 else "")
    basefolder = '/home2/qliu10/Transignaler_bp';
    print(basefolder, para_dict['save_file_pref'])

    '''
    para_dict['save_frequency'] = 10000;
    para_dict['train_epoch']  =1;

    para_dict['warmup_seq'] = 10000
    para_dict['save_frequency_s'] = para_dict['warmup_seq']//10;
    para_dict['save_frequency'] = para_dict['save_frequency_s']*5
    para_dict['save_frequency_s'] = para_dict['warmup_seq']//50;
    para_dict['save_frequency'] = para_dict['save_frequency_s']*5
    para_dict['train_epoch']  = 3;
    para_dict['cpu_devices'] = 7;
    para_dict['cuda_devices'] = 1;
    basefolder = '/data/data1/NanoporeRNA/Transignaler_bp'
    #finetuneTrainer( para_dict = para_dict, datafolders=basefolder+'/test_tfrecords' )
    '''
    # 337,837
    para_dict['warmup_seq'] = 10000
    para_dict['save_frequency'] = 10000 ; #para_dict['warmup_seq']
    para_dict['log_frequency']  =   100 ; #para_dict['warmup_seq']//10

    para_dict['train_epoch']  = 50
    #para_dict['lossfun'] = sys.argv[1]
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
    if len(sys.argv)>8:
        para_dict['ft_learning_rate'] = float(sys.argv[8])   
        para_dict['save_file_pref'] = para_dict['save_file_pref'] + ".ft"+str(int(float(sys.argv[8])*1e7+0.5))
    if len(sys.argv)>9:
        para_dict['saved_model'] = sys.argv[9]
    

    print(datetime.datetime.now())
    try:
        is_ft_t = (True if sys.argv[4] in [1, '1', 'T', 'True', 'true'] else False)
        datafs = [basefolder+'/tfrecord_ft.5M/Curlcake_m6a/', basefolder+'/tfrecord_ft.5M/Curlcake_non/']
        #finetuneTrainer( para_dict = para_dict, datafolders=datafs, is_ft=is_ft_t, ref_seq_file='EpiNano/Reference_sequences/cc.fasta',index_file="tfrecord.index", adam_learning_rate=int(sys.argv[5])/1e7, input_batch_size=int(sys.argv[6]), random_seed=3, not_use_chr=['cc6m_2244_T7_ecorv'] , downsampling_rate=(float(sys.argv[7])if len(sys.argv)>7 else None), saved_model=(sys.argv[9] if len(sys.argv)>9 else None))
        finetuneTrainer( para_dict = para_dict, datafolders=datafs, is_ft=is_ft_t, ref_seq_file='EpiNano/Reference_sequences/cc.fasta',index_file="tfrecord.index", random_seed=3, not_use_chr=['cc6m_2244_T7_ecorv'] )
    except Exception as e:
        print(traceback.format_exc())
        print(e)
        print(datetime.datetime.now())
    print(datetime.datetime.now())


