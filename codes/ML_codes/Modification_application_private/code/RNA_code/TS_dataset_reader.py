
import os,sys

import numpy as np
import random

import copy

import torch
import tensorflow as tf

sys.path.append( os.path.dirname(__file__) )

from TS_global import *
#import generate_tf_feature
import generate_tfrecord

class TS_Dataset(torch.utils.data.Dataset):
    def __init__(self , dataset_folder, index_file='ts_index.idx', isgzip=True, m_seed=1):
        self.dataset_folder = dataset_folder
        self.index_file = index_file;

        with open(self.dataset_folder+'/'+self.index_file, 'r') as f:
            self.tffiles = f.readlines()

        self.isgzip = isgzip;
        self.input_len = input_len
        self.pad_index = 0;
        self.m_seed = m_seed

        torch.manual_seed(self.m_seed)
        np.random.seed(self.m_seed)
        random.seed(self.m_seed)

        self.total_seq_count = 0;
        self.total_file_seq_count = 0;
        for _fl in self.tffiles:
            _fl_sp = _fl.strip().split()
            if len(_fl_sp)==1: continue;
            try:
                self.total_seq_count += int(_fl_sp[1])
                self.total_file_seq_count += 1;
            except:
                pass;

    def __len__(self):
        return len(self.tffiles)

    # https://github.com/pytorch/pytorch/issues/13246
    # use deepcopy in __getitem__ can solve this memory leaking problem.
    def __getitem__(self, idx):
        #this_file = self.dataset_folder+'/'+self.tffiles[idx].strip()
        this_file = self.tffiles[idx].strip().split()[0]
       
        tf_np_data = []
        tf_m_data = []
        tf_len_data = []

        #print("Process: open {}".format( this_file ))
        if self.isgzip:
            tf_data = tf.data.TFRecordDataset(this_file, compression_type='GZIP')
        else: tf_data = tf.data.TFRecordDataset(this_file)
        #tf_data.close();
        #print("Process: map {}".format( this_file ))
        tf_dataset = tf_data.map(generate_tfrecord.get_feature_from_tfrec)
        #print("Process: {}".format( this_file ))
        #procci = 0

        for tf_feat in tf_dataset.take(-1):
            tf_np_data.append( tf_feat[signal_features_in_tf].numpy())
            tf_m_data.append( tf_feat[mean_features_in_tf].numpy())
            #tf_len_data.append( tf_feat[len_features_in_tf].numpy())
            #if len( tf_np_data )>1: break;
            #if len( tf_np_data )<10:
            #    print( type( tf_np_data[-1] ), len(tf_np_data[-1]))
            #print( type(tf_np_data[-1]), len(tf_np_data[-1]))
            """
            procci += 1
            print( type(tf_np_data[-1] ), tf_np_data[-1].shape )
            random_center_pos = np.random.randint(1, input_len//2);
            for _testi in range ( len(tf_np_data[-1]) ):
                if _testi>100: break;
                outpt_s = ['{: 5d}'.format( _testi )];
                outpt_s.append( " | {:.3f}".format( tf_m_data[-1][_testi] ));
                for _idop in range( len(tf_np_data[-1][_testi]) ):
                    if tf_np_data[-1][_testi][_idop]<0.0000001: outpt_s.append('');
                    else: 
                        #print( type(tf_np_data[-1][_testi][_idop] ) , tf_np_data[-1][_testi][_idop] ) 
                        outpt_s.append("{:.1f}".format( tf_np_data[-1][_testi][_idop] ));
                print(';'.join(outpt_s))
            for _this_cp in range(random_center_pos, len(tf_np_data[-1]), input_len):
                _this_start = _this_cp - (input_len//2)
                _this_end = _this_cp + (input_len//2)+1

                meaning_start = dis_to_tail if _this_start>=0 else -_this_start + dis_to_tail
                meaning_end = input_len - dis_to_tail - (0 if _this_end<len(tf_np_data[-1]) else (_this_end-len(tf_np_data[-1])) ); # input_len - 2 if _this_end<len(tf_np_data[-1]) else (input_len-2-(_this_end-len(tf_np_data[-1])) )
                if _this_cp<=100: print(_this_start, _this_end, _this_cp, input_len, len(tf_np_data[-1]), input_len - dis_to_tail - (0 if _this_end<len(tf_np_data[-1]) else (_this_end-len(tf_np_data[-1])) ), meaning_start, meaning_end  )
                t_mask_pos = np.random.randint( meaning_start, meaning_end)

                arind = 0;
                while _this_start<0:
                    if _this_cp<=100: print( "s-{} {}".format( arind, _this_start ) )
                    _this_start +=  1;
                    arind += 1;
                while _this_start < _this_end:
                    if _this_start>=len(tf_np_data[-1]):
                        if _this_cp<=100: print( "s-{} {}".format( arind, _this_start ) )
                    else:
                        if _this_cp<=100:
                            outpt_s = ['s-{} {: 5d}'.format( arind,_this_start )];
                            outpt_s.append( " {} {:.3f}".format( "*" if (arind==t_mask_pos) else "|", tf_m_data[-1][_this_start] ));
                            for _idop in range( len(tf_np_data[-1][_this_start]) ):
                                if tf_np_data[-1][_this_start][_idop]<0.0000001: outpt_s.append('');
                                else:
                                    outpt_s.append("{:.1f}".format( tf_np_data[-1][_this_start][_idop] ));
                            print(';'.join(outpt_s))
                    _this_start +=  1;
                    arind += 1;

                if _this_cp<=100:
                    print("\t{} {}/ {} {} | {:.3f}".format(procci, t_mask_pos, _this_cp - (input_len//2), t_mask_pos+_this_cp - (input_len//2), tf_m_data[-1][ t_mask_pos+_this_cp - (input_len//2) ] ))
            """
        #print("load 1")

        return (tf_np_data, tf_m_data, tf_len_data);

if __name__=='__main__':
    cpuCount = os.cpu_count()

    # Print the number of
    # CPUs in the system
    print("Number of CPUs in the system:", cpuCount)
    basefolder = '/home2/qliu10/Transignaler_bp'; #'/data/data1/NanoporeRNA/Transignaler_bp'
    index_file = 'ts_record.idx'
    index_file = 'ts_record.s.idx'
    index_file = 'tfrecord.index'
    ds_folder = basefolder+'/test_tfrecords'; #'/data/data1/NanoporeRNA/tfrecords'
    
    ds_folder = basefolder+'/tfrecord.5M/na12878_rna_IVT/'
    index_file = 'tfrecord.index'

    signal_limits=m_md_norm_limits

    for train_epoch_i in range(3):
        test_ds = TS_Dataset(ds_folder, index_file, m_seed=train_epoch_i)

        test_datal = torch.utils.data.DataLoader(test_ds, batch_size=1, num_workers=1)

        test_get = next(iter(test_datal))
        print( type(test_get), len(test_get))
        for _i in range(3):
            print( len(test_get[_i]), type(test_get[_i]), (test_get[_i][0].shape if len(test_get[_i])>3 else "") )
        print('\n')
      
        if train_epoch_i==2:
            # t_ind = int( (_cs - signal_limits[0])/signal_limits[2] +0.5 )
            mea_mlist = []
            cal_mlist = []
            for seqid in range(len( test_get[1] ) ):
                #seqid = 0;
                testi = test_get[1][seqid].shape[1]//2;
                printlist = ["{:.3f}".format (test_get[1][seqid][0][testi] ) ]
                mea_mlist.append( test_get[1][seqid][0][testi] );
                cal_m = 0 
                for _vind in range( test_get[0][seqid].shape[2] ):
                    if test_get[0][seqid][0][ testi ][_vind]>1e-3:
                        printlist.append( "{}:{}".format( _vind, test_get[0][seqid][0][ testi ][_vind] ) )
                        cal_m += test_get[0][seqid][0][ testi ][_vind] * (_vind*signal_limits[2]+signal_limits[0])
                printlist.append( "{:.3f}".format (cal_m/100.0 ))
                cal_mlist.append( cal_m/100.0 )
                print(" ".join( printlist ))
            
            print('\n')
            dif_list = np.subtract(mea_mlist,cal_mlist)
            diflistop = []
            for _dv in dif_list:
                diflistop.append( "{:.3f}".format( _dv ) );
            print( len(mea_mlist), len(cal_mlist), np.square(np.subtract(mea_mlist,cal_mlist)).mean(), ', '.join(diflistop)  )
"""
<class 'list'> 3
543 <class 'list'> torch.Size([1, 118, 100])
543 <class 'list'> torch.Size([1, 118])
0 <class 'list'>
"""
