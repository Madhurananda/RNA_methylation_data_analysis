
m_md_norm_limits = [-5, 5, 0.1];
m_pa_norm_limits = [50, 150, 1]

dis_size = 100;

ValueCloseToZero = 1e-10/dis_size

input_len = 51;
input_len = 31;

#min_seq_len = 100;

dis_to_tail = 0;
dis_to_tail = 3;

signal_features_in_tf = 'dis_signals'
mean_features_in_tf = 'mean_signals'
len_features_in_tf = "len_signals"
fn_features_in_tf = "f5n"
bp_features_in_tf = 'bp'

TF_Withfn = False
TF_Withlen = False;
pos_features_in_tf = 'alignment_pos'

#mapped_map_info = 'map_info'
mapped_chrom_in_tf = 'mapped_chrom'
mapped_end_in_tf = 'mapped_end' 
mapped_start_in_tf = 'mapped_start' 
mapped_strand_in_tf = 'mapped_strand'


#                      ,  ,
tf_record_max_event = 1000000
tf_record_max_event =  500000


g_dropout_rate = 0.2
#g_dropout_rate = 0.1

#bp_to_int_dict = {'A':0, 'C':1, 'G':2, 'T':3, \
#                  'a':0, 'c':1, 'g':2, 't':3}

