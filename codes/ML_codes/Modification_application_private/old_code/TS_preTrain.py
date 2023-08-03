
import os,sys
import random
import time
import numpy as np
import copy
import math
import torch
import torch.multiprocessing

sys.path.append( os.path.dirname(__file__) )

import TS_model
import TS_dataset_reader

import TS_global

class train_feat_Batch:
    def __init__(self, data):
        transposed_data = list(zip(*data))
        self.feat = torch.stack(transposed_data[0], 0)
        self.pad = torch.stack(transposed_data[1], 0)
        self.masked = torch.stack(transposed_data[2], 0)
        self.maskLab = torch.stack(transposed_data[3], 0)

    def pin_memory(self):
        self.feat = self.feat.pin_memory()
        self.pad = self.pad.pin_memory()
        self.masked = self.masked.pin_memory()
        self.maskLab = self.maskLab.pin_memory()
        return self

def train_feat_Batch_wrapper(batch):
    return train_feat_Batch(batch)

def generate_train_feat(signal_full_info, p_batch_size, used_cpu_devices):
    feat_ts = []
    pad_inf = []
    mask_inf = []
    mask_label = []

    total_signal_event = 0
    signal_full_length = copy.deepcopy(signal_full_info[0])
    for _i_seq in range(len(signal_full_length)):
        total_signal_event += signal_full_length[_i_seq].size(1)
        if signal_full_length[_i_seq].size(1)<TS_global.min_seq_len: continue;
        random_center_pos = np.random.randint(1, TS_global.input_len//2);
        for _this_cp in range(random_center_pos, signal_full_length[_i_seq].size(1), TS_global.input_len):
            this_feat = []
            this_pad = []
            _this_start = _this_cp - (TS_global.input_len//2)
            _this_end = _this_cp + (TS_global.input_len//2)+1

            meaning_start = TS_global.dis_to_tail if _this_start>=0 else -_this_start + TS_global.dis_to_tail
            meaning_end = TS_global.input_len - TS_global.dis_to_tail - (0 if _this_end<signal_full_length[_i_seq].size(1) else (_this_end-signal_full_length[_i_seq].size(1)) ); #input_len - 2 if _this_end<signal_full_length[_i_seq].size(1) else (_this_end-signal_full_length[_i_seq].size(1)-2)
            t_mask_pos = np.random.randint( meaning_start, meaning_end)
            mask_inf.append( t_mask_pos )

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

            mask_label.append( this_feat[ t_mask_pos ]/100.0 + TS_global.ValueCloseToZero )
            #this_feat[ t_mask_pos ] = torch.zeros( TS_global.dis_size )
            random_substitute = np.random.randint(1, 100);
            if random_substitute<70:
                this_feat[ t_mask_pos ] = torch.zeros( TS_global.dis_size )
            elif random_substitute<85:
                this_feat[ t_mask_pos ] = torch.rand(TS_global.dis_size)
                this_feat[ t_mask_pos ] = this_feat[ t_mask_pos ]*100/torch.sum(this_feat[ t_mask_pos ] );
            feat_ts.append( torch.stack(this_feat) )
            pad_inf.append( this_pad )

    pad_inf = torch.tensor(pad_inf)
    mask_inf = torch.tensor(mask_inf)#.unsqueeze(1)
    feat_ts = torch.stack(feat_ts)
    mask_label = torch.stack (mask_label)

    train_dataset = torch.utils.data.TensorDataset(feat_ts, pad_inf, mask_inf, mask_label) ; #, id_info)
    train_loader = torch.utils.data.DataLoader(train_dataset, batch_size=p_batch_size, collate_fn=train_feat_Batch_wrapper, pin_memory=True, shuffle=True, num_workers=used_cpu_devices)

    return train_loader;


def get_adapt_learning_rate(init_lr, warmup_seq, num_train_seq, max_estimat_seq):
    decay_rate = 1;
    if num_train_seq <= warmup_seq:
        decay_rate = num_train_seq/warmup_seq
    else: # default learning rate is very small
        #if num_train_seq>max_estimat_seq*0.99: #-max_estimat_seq*0.01:
        if num_train_seq>(max_estimat_seq-warmup_seq)*0.99:
            #decay_rate = np.min([1- max_estimat_seq*0.01/(max_estimat_seq-warmup_seq), 
            #                 np.power( (max_estimat_seq*0.01/(warmup_seq)), -0.5 )]);
            #decay_rate = math.cos( math.pi/2*max_estimat_seq*0.99/(max_estimat_seq-warmup_seq) )
            decay_rate = math.cos( math.pi/2*(max_estimat_seq-warmup_seq)*0.99/(max_estimat_seq-warmup_seq) )
        else:
            #decay_rate = np.min([1- ( num_train_seq-warmup_seq )/(max_estimat_seq-warmup_seq),
            #                np.power( ((num_train_seq-warmup_seq )/(warmup_seq)), -0.5 )]);
            decay_rate = math.cos( ( math.pi/2*(num_train_seq-warmup_seq) )/(max_estimat_seq-warmup_seq) )
    return decay_rate*init_lr

def save_pretrained_model(mask_signal_model, adam_optmizer, device, para_dict, epoch, num_train_seq, loss_func, loss_func0, loss_func1, num_train_1epoch, max_estimat_seq):
    output_path = para_dict['save_folder']+'/'+para_dict['save_file_pref'] + "GPU.ep{}.{}.pt".format(epoch, num_train_seq)
    torch.save({
                "epoch": epoch,
                "num_train_seq": num_train_seq,
                "model_state_dict": mask_signal_model.state_dict(),
                "optimizer_state_dict": adam_optmizer.state_dict(),
                "loss_func": loss_func,
                "loss_func0": loss_func0, 
                "loss_func1": loss_func1,
                "num_train_1epoch": num_train_1epoch,
                "max_estimat_seq": max_estimat_seq
                }, 
                output_path)
    if 'cpumodel' in para_dict:
        output_path = para_dict['save_folder']+'/'+para_dict['save_file_pref'] + "CPU.ep{}.{}.pt".format(epoch, num_train_seq)
        torch.save({
                    "epoch": epoch,
                    "num_train_seq": num_train_seq,
                    "model_state_dict": mask_signal_model.cpu().state_dict(), 
                    "optimizer_state_dict": adam_optmizer.state_dict(),
                    "loss_func": loss_func,
                    "loss_func0": loss_func0,
                    "loss_func1": loss_func1,
                    "num_train_1epoch": num_train_1epoch,
                    "max_estimat_seq": max_estimat_seq
                    }, 
                    output_path)
        mask_signal_model.to(device)
        #adam_optmizer.to(device)
    print("EP:{}/{} Model Saved on: {}".format( epoch, num_train_seq, output_path));

# https://arxiv.org/pdf/1711.00489.pdf
# DON’T DECAY THE LEARNING RATE, INCREASE THE BATCH SIZE
# https://arxiv.org/pdf/1705.08741.pdf
# when using large batch sizes there is a persistent degradation in generalization performance - known as the "generalization gap" phenomenon. 
# 
#
#def preTrainer(para_dict, datafolder, index_file='tfrecord.index', input_batch_size: int=2048, cuda_devices: int=1, cpu_devices: int=5, adam_learning_rate: float =1e-4, adam_betas=(0.9,0.999), adam_weight_decay: float=0.01, random_seed: int=1,  train_epoch: int=10, warmup_seq: int=100000, gpu_saved_model=None):
def preTrainer(para_dict, datafolder, index_file='tfrecord.index', input_batch_size: int=2048, cuda_devices: int=1, cpu_devices: int=5, adam_learning_rate: float =1e-5, adam_betas=(0.9,0.999), adam_weight_decay: float=0.01, random_seed: int=1,  train_epoch: int=60, warmup_seq: int=1000000, gpu_saved_model=None):
    torch.manual_seed(random_seed)
    np.random.seed(random_seed)
    random.seed(random_seed)

    torch.multiprocessing.set_sharing_strategy('file_system')

    if not os.path.isdir(para_dict['save_folder']):
        os.system('mkdir -p {}'.format(para_dict['save_folder']))

    seqstep = para_dict['save_frequency']
    logstep = para_dict['log_frequency']
    if 'gpu_saved_model' in para_dict:
        gpu_saved_model = para_dict['gpu_saved_model']
    if 'warmup_seq' in para_dict:    
        if para_dict['warmup_seq']<1000: warmup_seq    = para_dict['warmup_seq']
        else: print('warmup_seq is too small: {}. not used.'.format(para_dict['warmup_seq']))
    if 'train_epoch' in para_dict:   
        if para_dict['train_epoch']>0: train_epoch   = para_dict['train_epoch']
        else: print('train_epoch is too small: {}. not used.'.format(para_dict['train_epoch']))
    if 'random_seed' in para_dict:   
        random_seed   = para_dict['random_seed']
    if 'cpu_devices' in para_dict:   
        if para_dict['cpu_devices']>=0: cpu_devices   = para_dict['cpu_devices']
        else: print('cpu_devicesis too small: {}. not used.'.format(para_dict['cpu_devices']))
    if 'cuda_devices' in para_dict:  
        if para_dict['cuda_devices']>0: cuda_devices  = para_dict['cuda_devices']
        else:  print('cuda_devices too small: {}. not used.'.format(para_dict['cuda_devices']))
    if 'input_batch_size' in para_dict: 
        if para_dict['input_batch_size']>0: input_batch_size = para_dict['input_batch_size']
        else: print('input_batch_sizeis too small: {}. not used.'.format(para_dict['input_batch_size']))
    if 'index_file' in para_dict:    
        index_file    = para_dict['index_file']

    print("Key learn info input_batch_size={} adam_learning_rate={} warmup_seq={:,} seqstep={:,} logstep={:,}".format(input_batch_size, adam_learning_rate, warmup_seq, seqstep, logstep))

    mask_signal_model = TS_model.SignalBertPred( ispretrain = True );
    device = torch.device("cuda:0" if torch.cuda.is_available() else "CPU")
    mask_signal_model = mask_signal_model.to(device);
    para_sum = 0;
    for l_name, l_para in mask_signal_model.named_parameters():
        #print(l_name, l_para.data.shape, l_para.requires_grad,  type(l_name) )
        if l_name.split('.')[0] not in ['pred_linear', 'pred_linear1', 'pred_linear2']:
            para_sum += l_para.nelement(); #l_para.data.size(0);
    print("Total parameters: {}".format( para_sum ))   

    cpuCount = os.cpu_count()
    used_cpu_devices = cpu_devices
    if cpuCount <= used_cpu_devices : used_cpu_devices = cpuCount
    print("Number of CPUs in the system: {}. Used CPUs: {}".format(cpuCount, used_cpu_devices ) )

    if torch.cuda.device_count() > 1 and cuda_devices>0:
        if cuda_devices>=torch.cuda.device_count(): used_cuda_devices = torch.cuda.device_count()
        else: used_cuda_devices = cuda_devices
        print("Number of GPUs in the system: {}. {} GPUs for training".format( torch.cuda.device_count(), used_cuda_devices ))
        #if used_cuda_devices>1: 
        mask_signal_model = torch.nn.DataParallel( mask_signal_model, device_ids=list(range( used_cuda_devices )) );

    adam_optmizer = torch.optim.Adam( mask_signal_model.parameters(), lr = adam_learning_rate, betas=adam_betas, weight_decay=adam_weight_decay )
    loss_func  = torch.nn.KLDivLoss(reduction="batchmean");
    loss_func0 = torch.nn.KLDivLoss(reduction="batchmean");
    loss_func1 = torch.nn.KLDivLoss(reduction="batchmean");
    #llog_func = torch.nn.KLDivLoss(reduction="batchmean", log_target=True)
    #loss_cels = torch.nn.CrossEntropyLoss()
   
    max_estimat_seq = 0;
    num_train_seq = 0;
    save_model_f = 0;
    num_train_i = 0
    first_n_warmup = 3;
    train_epoch_i = 0
    num_train_1epoch = 0;

    if not gpu_saved_model==None:
        t_checkpoint = torch.load( gpu_saved_model); #, map_location=device );
        mask_signal_model.load_state_dict(t_checkpoint['model_state_dict'])
        adam_optmizer.load_state_dict(t_checkpoint['optimizer_state_dict'])
        num_train_seq = t_checkpoint['num_train_seq']
        save_epoch_i = t_checkpoint['epoch']
        train_epoch_i = t_checkpoint['epoch']-1
        loss_func = t_checkpoint['loss_func']
        loss_func0 = t_checkpoint['loss_func0']
        loss_func1 = t_checkpoint['loss_func1']
        num_train_1epoch = t_checkpoint['num_train_1epoch']
        max_estimat_seq = t_checkpoint['max_estimat_seq']

        total_max_estimat_seq = max_estimat_seq*(1 if train_epoch<2 else train_epoch//2)

        load_num_train_seq = 0;
        for _prev_epoch_i in range(0, train_epoch_i):
            load_num_train_seq += num_train_1epoch

        """fs_dataset = TS_dataset_reader.TS_Dataset(datafolder, index_file, train_epoch_i+1)
        fs_dataloader = torch.utils.data.DataLoader(fs_dataset, batch_size=1, num_workers=used_cpu_devices)
        for batch_fidx, batch_data in enumerate( fs_dataloader ):
            if load_num_train_seq+len(batch_data[0])>=num_train_seq:
                break;
            load_num_train_seq += len(batch_data[0])"""

    #avg_los = []
    sum_los = 0;
    previous_num_train_i = num_train_i
    start_time = time.time()
    #for train_epoch_i in range(train_epoch):
    while (train_epoch_i < train_epoch) or (train_epoch_i>=train_epoch and  new_lr>1e-8) :
        train_epoch_i += 1; # epoch_i starts from 1 to train_epoch
        
        #print(train_epoch_i, save_epoch_i, train_epoch)
        #if gpu_saved_model==None or train_epoch_i>save_epoch_i:
        fs_dataset = TS_dataset_reader.TS_Dataset(datafolder, index_file, m_seed=train_epoch_i)
        if train_epoch_i==1: 
            if fs_dataset.total_file_seq_count>=len(fs_dataset)-10 and fs_dataset.total_file_seq_count>1 :
                max_estimat_seq = fs_dataset.total_seq_count
            else:
                max_estimat_seq = len(fs_dataset)*(800 if TS_global.tf_record_max_event ==1000000 else 400)
            total_max_estimat_seq = max_estimat_seq*(1 if train_epoch<2 else train_epoch//2)
            print(type(fs_dataset), len(fs_dataset), fs_dataset.total_file_seq_count, fs_dataset.total_seq_count )
            sys.stdout.flush()
        fs_dataloader = torch.utils.data.DataLoader(fs_dataset, batch_size=1, num_workers=used_cpu_devices)
        
        for batch_fidx, batch_data in enumerate( fs_dataloader ):
            t_seq_num = len(batch_data[0])
            if (not gpu_saved_model==None) and load_num_train_seq+t_seq_num<=num_train_seq:
                load_num_train_seq +=t_seq_num
                continue;
            #print( batch_fidx, len(batch_data), type(batch_data), "Detail0:",  len(batch_data[0]), len(batch_data[0][0]), len(batch_data[0][0][0]), len(batch_data[0][0][0][0]), type(batch_data[0][0][0][0][0]) , type(batch_data[0]), batch_data[0][0].shape, "1/2", len(batch_data[1]), type(batch_data[1]), len(batch_data[2]), type(batch_data[2]))
            #for _testi in range(10):
            #        print("\t",  batch_data[0][_testi].shape);
            save_model_f = 0; 
            # seqstep=10,000, train_epoch=10, warmup_seq=100,000
            # seqstep=20,000, train_epoch=3,  warmup_seq=100,000 ; 20,000,000 = 200
            #if (train_epoch_i==1 and ((num_train_seq+t_seq_num)//(logstep if num_train_seq<warmup_seq*first_n_warmup else seqstep)) > (num_train_seq//(logstep if num_train_seq<warmup_seq*first_n_warmup else seqstep)) ) or \
            #        ((num_train_seq+t_seq_num)//(seqstep*10)) > (num_train_seq//(10*seqstep)):
            if ((num_train_seq+t_seq_num)//logstep) > (num_train_seq//logstep):
                save_model_f = 1;
            if ((num_train_seq+t_seq_num)//seqstep) > (num_train_seq//seqstep):
                save_model_f = 2
            num_train_seq += t_seq_num
            if train_epoch_i==1:
                num_train_1epoch += t_seq_num
            new_lr = get_adapt_learning_rate(adam_learning_rate, warmup_seq, num_train_seq, total_max_estimat_seq);
            for param_group in adam_optmizer.param_groups:
                param_group['lr'] = new_lr

            train_loader = generate_train_feat( batch_data, input_batch_size, used_cpu_devices)
            del batch_data
            for batch_ndx, train_batch in enumerate(train_loader):
                #print( train_batch.feat.shape, train_batch.masked.shape, train_batch.pad.shape, train_batch.feat.shape)
                if train_batch.feat.size(0) < input_batch_size/2: continue;
                num_train_i += 1

                tb_maskpos = train_batch.masked.to(device)
                tb_pad = train_batch.pad.to(device)
                tb_feat = train_batch.feat.to(device)
                pred_mask, pred_logmask = mask_signal_model.forward(tb_feat, tb_pad, device, tb_maskpos)
                if torch.isnan(pred_logmask).any() or torch.isinf(pred_logmask).any():
                    print("Warning!!!! nan or inf value in log-pred {}".format(pred_logmask))
                    continue;
                aggrat_pred_info = [ torch.zeros(pred_mask.size(0), pred_mask.size(1)//2), torch.zeros(pred_mask.size(0), pred_mask.size(1)//2+1)]
                aggrat_labl_info = [ torch.zeros(pred_mask.size(0), pred_mask.size(1)//2), torch.zeros(pred_mask.size(0), pred_mask.size(1)//2+1)]
                for _c_ind in range(0, pred_mask.size(1), 2):
                    aggrat_pred_info[0][:,_c_ind//2] = (pred_mask[:, _c_ind] + pred_mask[:, _c_ind+1]).log(); 
                    aggrat_labl_info[0][:,_c_ind//2] = (train_batch.maskLab[:, _c_ind] + train_batch.maskLab[:, _c_ind+1]); #.log();
                aggrat_pred_info[1][:,0] = pred_mask[:,0].log();
                aggrat_labl_info[1][:,0] =train_batch.maskLab[:, 0]; #.log();
                for _c_ind in range(1, pred_mask.size(1), 2):
                    aggrat_pred_info[1][:,_c_ind//2+1] = ((pred_mask[:, _c_ind] + pred_mask[:, _c_ind+1]) if _c_ind+1 <pred_mask.size(1) else pred_mask[:, _c_ind]).log();
                    aggrat_labl_info[1][:,_c_ind//2+1] =  (train_batch.maskLab[:, _c_ind] + train_batch.maskLab[:, _c_ind+1]) if _c_ind+1 <pred_mask.size(1) else train_batch.maskLab[:, _c_ind];

                if torch.isnan(aggrat_pred_info[0]).any() or torch.isnan(aggrat_pred_info[1]).any() or \
                   torch.isinf(aggrat_pred_info[0]).any() or torch.isinf(aggrat_pred_info[1]).any() :
                    if torch.isnan(aggrat_pred_info[0]).any() or torch.isinf(aggrat_pred_info[0]).any():
                       print("Warning!!!! nan or inf value in log-pred {}".format( aggrat_pred_info[0] ))
                    if torch.isnan(aggrat_pred_info[1]).any() or torch.isinf(aggrat_pred_info[1]).any():
                       print("Warning!!!! nan or inf value in log-pred {}".format( aggrat_pred_info[1] ))
                    continue;

                loss = loss_func(pred_logmask, train_batch.maskLab.to(device)) + \
                       loss_func0(aggrat_pred_info[0].to(device), aggrat_labl_info[0].to(device) ) + \
                       loss_func1(aggrat_pred_info[1].to(device), aggrat_labl_info[1].to(device) )
    
                adam_optmizer.zero_grad(); 
                loss.backward();
                adam_optmizer.step();
                #avg_los.append(loss.item());
                sum_los += loss.item()

                if num_train_i<10:
                    print("This loss: {:.3f} for ep{}.{:,}. Consuming time: {:.1f}. Learning_rate={:.5e} avgloss={:.3f} traini={:,}".format(loss.item(), train_epoch_i, num_train_seq, time.time()-start_time, new_lr, sum_los/(num_train_i-previous_num_train_i), num_train_i ))

            # end for batch_ndx, train_batch in enumerate(train_loader):
            """if save_model_f==0 and (num_train_seq<warmup_seq//2):
                #print("This loss: {:.3f} for ep{}.{:,}. Consuming time: {:.1f}. Learning_rate={:.5e} avgloss={:.3f}/{:.3f}".format(loss.item(), train_epoch_i, num_train_seq, time.time()-start_time, new_lr, np.mean(avg_los), np.std(avg_los) ))
                print("This loss: {:.3f} for ep{}.{:,}. Consuming time: {:.1f}. Learning_rate={:.5e} avgloss={:.3f} traini={:,}".format(loss.item(), train_epoch_i, num_train_seq, time.time()-start_time, new_lr, sum_los/(num_train_i-previous_num_train_i), num_train_i ))
                start_time = time.time()
                sum_los =0;
                previous_num_train_i = num_train_i
            if save_model_f>0:"""
            if save_model_f>0 or (num_train_seq<warmup_seq//2):
                #print("This loss: {:.3f} for ep{}.{:,}. Consuming time: {:.1f}. Learning_rate={:.5e} avgloss={:.3f}/{:.3f}".format(loss.item(), train_epoch_i, num_train_seq, time.time()-start_time, new_lr, np.mean(avg_los), np.std(avg_los) )) #, loss_func( mask_signal_model.forward(tb_feat, tb_pad, device, tb_maskpos)[0] , train_batch.maskLab.to(device)).item() ))
                print("This loss: {:.3f} for ep{}.{:,}. Consuming time: {:.1f}. Learning_rate={:.5e} avgloss={:.3f} traini={:,}".format(loss.item(), train_epoch_i, num_train_seq, time.time()-start_time, new_lr, sum_los/(num_train_i-previous_num_train_i), num_train_i ))
                start_time = time.time()
                #avg_los = []
                sum_los =0;
                previous_num_train_i = num_train_i
                if save_model_f>1:
                    save_pretrained_model(mask_signal_model, adam_optmizer, device, para_dict, train_epoch_i, num_train_seq, loss_func, loss_func0, loss_func1, num_train_1epoch, max_estimat_seq)
                sys.stdout.flush()
        #end for batch_fidx, batch_data in enumerate( fs_dataloader ):
    save_pretrained_model(mask_signal_model, adam_optmizer, device, para_dict, train_epoch_i, num_train_seq, loss_func, loss_func0, loss_func1, num_train_1epoch, max_estimat_seq)

if __name__=='__main__':
    para_dict = {}
    para_dict['save_folder'] = './testModel';
    para_dict['save_file_pref'] = 'test.'
    basefolder = '/data/data1/NanoporeRNA/Transignaler_bp'
    """
    #                                 ,
    para_dict['log_frequency'] = 20000;
    para_dict['save_frequency'] = para_dict['log_frequency']*20
    #para_dict['save_frequency'] = 10
    para_dict['train_epoch']  = 3;
    para_dict['cpu_devices'] = 5;
    para_dict['cuda_devices'] = 2
    #para_dict['gpu_saved_model'] = './testModel/test.GPU.ep1.279.pt'

    #preTrainer( para_dict = para_dict, datafolder=basefolder+'/test_tfrecords', index_file="tfrecord_test.index" )
    #preTrainer( para_dict = para_dict, datafolder=basefolder+'/tfrecord1M', index_file="tfrecord_test.index" )
    #preTrainer( para_dict = para_dict, datafolder=basefolder+'/tfrecord1M', index_file="tfrecord.index" )
    preTrainer( para_dict = para_dict, datafolder=basefolder+'/tfrecord.5M', index_file="tfrecord.index" )
    """

    para_dict['log_frequency'] = 1000000//10;
    para_dict['save_frequency'] = para_dict['log_frequency']*5
    para_dict['warmup_seq'] = 1000000
    #para_dict['save_frequency'] = 10
    para_dict['train_epoch']  = 3;
    para_dict['cpu_devices'] = 7;
    para_dict['cuda_devices'] = 1;
    #para_dict['gpu_saved_model'] = './testModel/test.GPU.ep1.279.pt'
    
    #                           ,  ,
    para_dict['warmup_seq'] = 10000000
    #para_dict['log_frequency'] = para_dict['warmup_seq']//10;
    #para_dict['save_frequency'] = para_dict['log_frequency']*5
    #para_dict['log_frequency'] = para_dict['warmup_seq']//50;
    #para_dict['save_frequency'] = para_dict['log_frequency']*5
    para_dict['log_frequency'] = para_dict['warmup_seq']//200
    para_dict['save_frequency'] = para_dict['warmup_seq']/20
    para_dict['train_epoch']  = 20;
    para_dict['cpu_devices'] = 10;
    para_dict['cuda_devices'] = 1;

    #preTrainer( para_dict = para_dict, datafolder=basefolder+'/test_tfrecords', index_file="tfrecord_test.index" )
    #preTrainer( para_dict = para_dict, datafolder=basefolder+'/tfrecord1M', index_file="tfrecord_test.index" )
    #preTrainer( para_dict = para_dict, datafolder=basefolder+'/tfrecord1M', index_file="tfrecord.index" )
    #preTrainer( para_dict = para_dict, datafolder=basefolder+'/tfrecord.5M', index_file="tfrecord.index" )
    #preTrainer( para_dict = para_dict, datafolder=basefolder+'/tfrecord.5M', index_file="tfrecord.test.index" )
    print(time.time())
    try:
        preTrainer( para_dict = para_dict, datafolder=basefolder+'/tfrecord.5M', index_file="tfrecord.index", adam_learning_rate=1e-5 )
    except  Exception as e:
        print(e)
        print(time.time())


