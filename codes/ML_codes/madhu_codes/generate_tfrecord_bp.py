
import os,sys
import numpy as np
import h5py
import glob
import time

import tensorflow as tf
from tqdm import tqdm
import read_fasta
import read_tombo
import read_bam
import _generate_feature_

from TS_global import *


def create_tfrec(m_f, mean_list, len_list, f5n, bp_list):
    m_f = np.array(m_f); # m_f = m_f.reshape((-1, dis_size));
    #print(m_f.shape)
    #print(m_f[5,:])
    #global test_tflist
    #test_tflist.append(m_f[5,:])

    mean_array = np.array(mean_list)
    len_array = np.array(len_list)

    bp_array = np.array( bp_list )

    feature = {signal_features_in_tf : _generate_feature_.float_feature_list( m_f.flatten() ),
                mean_features_in_tf: _generate_feature_.float_feature_list(mean_array), 
                len_features_in_tf: _generate_feature_.float_feature_list(len_array),
                fn_features_in_tf: _generate_feature_.string_features(f5n),
                bp_features_in_tf: _generate_feature_.int64_feature_list(bp_array)
                }

    return tf.train.Example(features=tf.train.Features(feature=feature))

def get_feature_from_tfrec( tf_rec ): #, withfn=False, withlen=False ):
    feature_list = { signal_features_in_tf: tf.io.VarLenFeature(tf.float32), 
                    mean_features_in_tf: tf.io.VarLenFeature(tf.float32),
                    len_features_in_tf: tf.io.VarLenFeature(tf.float32),
                    fn_features_in_tf: tf.io.FixedLenFeature([], tf.string),
                    bp_features_in_tf: tf.io.VarLenFeature(tf.int64)
                    }
    m_feat = tf.io.parse_single_example(tf_rec, feature_list)
    m_feat[signal_features_in_tf] = tf.sparse.to_dense(m_feat[signal_features_in_tf])
    
    m_feat[mean_features_in_tf] = tf.sparse.to_dense(m_feat[mean_features_in_tf])
    m_feat[len_features_in_tf] = tf.sparse.to_dense(m_feat[len_features_in_tf])
    m_feat[bp_features_in_tf] = tf.sparse.to_dense(m_feat[bp_features_in_tf])

    m_feat[signal_features_in_tf] = tf.reshape(m_feat[signal_features_in_tf], [-1, dis_size])

    if not TF_Withfn:
        del m_feat[fn_features_in_tf]
    if not TF_Withlen:
        del m_feat[len_features_in_tf]

    return m_feat

def feat_generator(f5files_Q, tffiles_Q, RNA_DNA, tfrecords_dir, tf_t_record_max_event, f5_t_batch_idx, para_dict):
    f5r = read_tombo.Fast5Read()
    seqr = read_tombo.SeqReverse()

    tfrec_file_prex = "signal.tfrec"
    tfroptions = tf.io.TFRecordOptions(compression_type = 'GZIP')
    round_event_size = 0;
    round_times = 0;
    round_f5fnum = 0
    tf_writer = None; #tf.io.TFRecordWriter(tfrecords_dir+tfrec_file_prex++"_"+str(round_times)+".gz", options=tfroptions)
    previous_tf_fn = "";
    process_start_time = time.time();
    while not f5files_Q.empty():
        try:
            p_f5files, p_bath_idx = f5files_Q.get(block=False);
        except:
            break;

        if (p_bath_idx+1)%20==0:
            print("Process {} among {} batch with batch size {}. Processing time: {:.1f}.".format( p_bath_idx+1, f5_t_batch_idx, len(p_f5files), time.time()-process_start_time ))
            process_start_time = time.time();
            sys.stdout.flush()

        if tf_writer==None:
            tf_writer = tf.io.TFRecordWriter(tfrecords_dir+tfrec_file_prex+"_"+str(p_bath_idx)+"_"+str(round_times)+".gz", options=tfroptions)
            previous_tf_fn = tfrecords_dir+tfrec_file_prex+"_"+str(p_bath_idx)+"_"+str(round_times)+".gz"
        for f5f in (p_f5files):
        # for f5f in tqdm(p_f5files):
            if round_event_size >= tf_t_record_max_event: #tf_record_max_event:
                #tffiles_Q.put(tfrecords_dir+tfrec_file_prex+"_"+str(p_bath_idx)+"_"+str(round_times)+".gz" )
                tffiles_Q.put( (previous_tf_fn, round_f5fnum) )
                tf_writer.close()
                round_times += 1
                round_f5fnum = 0
                round_event_size = 0;
                tf_writer = tf.io.TFRecordWriter(tfrecords_dir+tfrec_file_prex+"_"+str(p_bath_idx)+"_"+str(round_times)+".gz", options=tfroptions)
                previous_tf_fn = tfrecords_dir+tfrec_file_prex+"_"+str(p_bath_idx)+"_"+str(round_times)+".gz"
            _mf, fsize, mean_list, len_list,  aligninfo, bp_list = _generate_feature_.generate_f( f5f, f5r, seqr, para_dict, RNA_DNA = RNA_DNA)
            if not fsize>para_dict['length_thr']: #min_seq_len: 
                continue;
            #if fsize>para_dict['length_thr']: #min_seq_len:
            tf_writer.write(create_tfrec(_mf, mean_list, len_list, f5f.split('/')[-1], bp_list ).SerializeToString() );
            round_f5fnum += 1;
            round_event_size += fsize
    if round_event_size>0:
        #tffiles_Q.put(tfrecords_dir+tfrec_file_prex+"_"+str(p_bath_idx)+"_"+str(round_times)+".gz" )
        tffiles_Q.put( (previous_tf_fn, round_f5fnum) )
    if not tf_writer==None: tf_writer.close()


def generate_f_parallel(f5filles, para_dict, numprocess=5, batch_file_size=4000, RNA_DNA='RNA', tfrecords_dir='test_tfrecordsbp/', tf_t_record_max_event=tf_record_max_event): #length_thr=100):
    import multiprocessing
    import time

    process_manager = multiprocessing.Manager();
    f5files_Q = process_manager.Queue();
    tffiles_Q = process_manager.Queue();

    if numprocess<1: numprocess = 1

    print("Total {} fast5 files".format( len(f5filles) ))
    sys.stdout.flush()

    f5batch = []
    f5_c_batch_idx = 0;
    for f5f in f5files:
        f5batch.append( f5f )
        if len(f5batch)==batch_file_size:
            f5files_Q.put((f5batch, f5_c_batch_idx))
            f5batch = [];
            f5_c_batch_idx += 1
    if len(f5batch)>0:
        f5files_Q.put((f5batch, f5_c_batch_idx))
        f5batch = [];
        f5_c_batch_idx += 1

    share_param = (f5files_Q, tffiles_Q, RNA_DNA, tfrecords_dir, tf_t_record_max_event, f5_c_batch_idx, para_dict)
    process_handlers = []
    for id_p in range(numprocess):
        c_proc = multiprocessing.Process(target=feat_generator, args=share_param);
        c_proc.start();
        process_handlers.append( c_proc )
    
    tfind_writer = open(tfrecords_dir+'/tfrecord.index', 'w')
    while any(c_p.is_alive() for c_p in process_handlers):
        try:
            tffn, round_f5fnum = tffiles_Q.get(block=False)
            tfind_writer.write(tffn+" "+str(round_f5fnum)+'\n')
        except:
            time.sleep(1)
            continue;
    while True:
        try:
            tffn, round_f5fnum = tffiles_Q.get(block=False)
            tfind_writer.write(tffn+" "+str(round_f5fnum)+'\n')
        except:
            break;
    tfind_writer.close()


if __name__ == '__main__':
    if len(sys.argv)<3: #    0  1  2  3  4  5  6  7  8  9                   0        1           2                3           4          5         6          7          8             9
        print("Usage: python {} {} {} {} {} {} {} {} {} {}".format( sys.argv[0], "f5folder", "output_folder", "#thread", "single_mult", "path", "bamfile", "len_thr", "map_qual_thr", "RD"))
        sys.exit(1)
    f5folder = sys.argv[1]
    
    # if len(sys.argv)>4 and int(sys.argv[4])>1:
    #     if len(sys.argv)>5 and len(sys.argv[5])>1:
    #         #f5files = glob.glob(os.path.join(f5folder, sys.argv[5]+"/*/*.fast5"));
    #         #print( os.path.join(f5folder, sys.argv[5], "*/*.fast5") )
    #         midpath = sys.argv[5]
    #         while len(midpath)>0:
    #             if midpath[0] in ['/', '\\']:
    #                 midpath = midpath[1:]
    #             else: break;
    #         f5files = glob.glob(os.path.join(f5folder, midpath, "*/*.fast5"));
    #     else:
    #         f5files = glob.glob(os.path.join(f5folder, "*/*.fast5"));
    # else:
    #     f5files = glob.glob(os.path.join(f5folder, "*.fast5"));
    
    f5files = [os.path.join(root, name)
             for root, dirs, files in os.walk(f5folder)
             for name in files
             if name.endswith((".fast5"))]

    tf_t_record_max_event = tf_record_max_event
    #if len(sys.argv)>6 and int(sys.argv[6])>tf_record_max_event/100:
    #    tf_t_record_max_event = int(sys.argv[6])
    tfrecords_dir = 'test_tfrecordsbp/'
    tfrecords_dir = sys.argv[2]
    if not os.path.isdir( tfrecords_dir ):
        os.system("mkdir -p {}".format(tfrecords_dir))

    numthread = int(sys.argv[3])
    # Number of CPUs in the system: 20. Used CPUs: 5
    # Number of GPUs in the system: 2. 1 GPUs for training

    ########### Do not comment this one off, I do not know why this one is here .... ##########
    # generate_tfrecord( f5files, tfrecords_dir )
    #######################################################################################

    para_dict = {}
    para_dict['length_thr'] = 100
    para_dict['length_thr'] = int(sys.argv[7])
    para_dict['map_quality'] = 20
    para_dict['map_quality'] = int(sys.argv[8])

    para_dict['bam_info'] = read_bam.read_bam( sys.argv[6], para_dict['map_quality'], para_dict['length_thr'] );

    generate_f_parallel( f5files, para_dict, numthread, tfrecords_dir=tfrecords_dir, tf_t_record_max_event=tf_t_record_max_event, RNA_DNA=sys.argv[9])




