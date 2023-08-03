
import os,sys
import numpy as np
import h5py
import glob

import read_fasta
import read_tombo

import tensorflow as tf


from TS_global import *

# https://keras.io/examples/keras_recipes/creating_tfrecords/
def int64_feature(value):
    """Returns an int64_list from a bool / enum / int / uint."""
    return tf.train.Feature(int64_list=tf.train.Int64List(value=[value]))
def int64_feature_list(value):
    """Returns an int64_list from a bool / enum / int / uint."""
    return tf.train.Feature(int64_list=tf.train.Int64List(value=value))
def float_feature_list(value):
    """Returns a list of float_list from a float / double."""
    return tf.train.Feature(float_list=tf.train.FloatList(value=value))

def string_features(value):
    return tf.train.Feature(bytes_list=tf.train.BytesList(value=[value.encode()]) )


def discritize(m_signal_seg, signal_limits=m_md_norm_limits):
    mdis_ = np.zeros(dis_size)

    signal_len = len(m_signal_seg)
    mod_ind = {}
    for _cs in m_signal_seg:
        t_ind = int( (_cs - signal_limits[0])/signal_limits[2] +0.5 )
        if t_ind<0: t_ind = 0;
        if t_ind>=dis_size: t_ind = dis_size-1
        mdis_[ t_ind ] += 1;
        mod_ind[ t_ind ] = True;
    for _mind in mod_ind:
        mdis_[_mind] = (mdis_[_mind]*100)//signal_len

    return mdis_

def generate_f(f5f, f5r, seqr, para_dict, signal_limits=m_md_norm_limits, RNA_DNA='RNA'):
    m_f = []
    mean_list = []
    len_list = []
    bp_list = []
    h5reader = h5py.File(f5f, 'r')
    #largdif = 0;
    tombo_sum = f5r.read_tombo_info( h5reader )
    if not f5r.tombo_no: #return (m_f, 0)              0                      1                   2              3              4               5
        tombo_align_summary = tombo_sum[0] ; # ["clipped_bases_end", "clipped_bases_start", "mapped_chrom", "mapped_end", "mapped_start", "mapped_strand"]
        tombo_scale_summary = tombo_sum[1] ; # ["lower_lim", "scale", "shift", "upper_lim"]

        fq_4 = f5r.read_fq( h5reader )
        good_align = False;
        read_id = fq_4[0].split()[0]
        if read_id[0] in ['@', '>']:
            read_id = read_id[1:]
        if read_id in para_dict['bam_info']:
            for c_align_inf in para_dict['bam_info'][read_id]:
                if tombo_align_summary[2]==c_align_inf[0] and ( (tombo_align_summary[5]=='+' and (not c_align_inf[2])) or (tombo_align_summary[5]=='-' and c_align_inf[2]) ):
                    good_align = True;
                    break;
        #else:
        if not good_align:
            print("Filter reads: {} {} {}".format( tombo_align_summary[2:], para_dict['bam_info'][read_id] if read_id in para_dict['bam_info'] else "Not-in", read_id ))

        if good_align:
            tombo_events = f5r.read_tombo_event( h5reader)
            _, normal_signals = f5r.read_signals(h5reader, tombo_scale_summary[1:3])
            if RNA_DNA=='RNA':
                normal_signals = normal_signals[::-1]
            for _c_e in range(len( tombo_events )):
                m_f.append(discritize(normal_signals[ tombo_events[_c_e][2]:(tombo_events[_c_e][2]+tombo_events[_c_e][3]) ], signal_limits))
                mean_list.append( np.mean( normal_signals[ tombo_events[_c_e][2]:(tombo_events[_c_e][2]+tombo_events[_c_e][3]) ] ) )
                len_list.append( tombo_events[_c_e][3] )
                #bp_list.append( bp_to_int_dict[ tombo_events[_c_e][4].decode('UTF-8') ] ) # add on Sep 21, 2022
                ## change on Feb 21, 2023
                bp_list.append( seqr.get_base_ind ( tombo_events[_c_e][4].decode('UTF-8') , tombo_align_summary[5] if not f5r.tombo_no else None, RNA_DNA ) )
                #if abs( tombo_events[_c_e][0] - np.mean( normal_signals[ tombo_events[_c_e][2]:(tombo_events[_c_e][2]+tombo_events[_c_e][3]) ]) )>0.005:
                #    largdif += 1;
                #    print("\t{} {:.3f} {:.3f} {:.3f}".format( _c_e, tombo_events[_c_e][0], np.mean( normal_signals[ tombo_events[_c_e][2]:(tombo_events[_c_e][2]+tombo_events[_c_e][3]) ]), tombo_events[_c_e][0]-np.mean( normal_signals[ tombo_events[_c_e][2]:(tombo_events[_c_e][2]+tombo_events[_c_e][3]) ])  ))
            #if largdif>0:
            #    print("LargeDif: {} {}/{}".format( f5f, largdif, len(tombo_events) ))
    h5reader.close()
    return (m_f, len(m_f), mean_list, len_list, (tombo_align_summary[2:] if not f5r.tombo_no else None), bp_list ); #torch.stack(m_f)

if __name__ == '__main__':
    f5folder = sys.argv[1]

    f5files = glob.glob(os.path.join(f5folder, "*.fast5"));




