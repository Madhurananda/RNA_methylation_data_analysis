
import os,sys
import numpy as np
import h5py
import glob

import read_fasta

class Fast5Read:
    def __init__(self):
        self.Analysis_path = '/Analyses'
        self.Basecall_path = 'Basecall_1D_000'
        self.template_path = 'BaseCalled_template'
        self.fastq_path = 'Fastq'
        self.move_path = 'Move'
        self.signal_path = '/Raw/Reads'
        self.tombo_path = 'RawGenomeCorrected_000'
        self.event_t_path = 'Events'
        self.align_t_path = 'Alignment'
        self.signal_entry = 'Signal'
        self.channel_path = '/UniqueGlobalKey/channel_id'

        self.signal_full_path = '/'.join([self.signal_path])
        self.fq_full_path = '/'.join([self.Analysis_path, self.Basecall_path, self.template_path, self.fastq_path])
        self.tombo_event_full_path = '/'.join([self.Analysis_path, self.tombo_path, self.template_path, self.event_t_path])
        self.tombo_align_full_path = '/'.join([self.Analysis_path, self.tombo_path, self.template_path, self.align_t_path])
        self.tombo_base_full_path = '/'.join([self.Analysis_path, self.tombo_path, self.template_path])
        self.tombo_align_attr = ["clipped_bases_end", "clipped_bases_start", "mapped_chrom", "mapped_end", "mapped_start", "mapped_strand"]

        self.channel_keyKeys = ['digitisation', 'offset', 'range', 'sampling_rate', 'channel_number']

        self.tombo_signal_norm_keys = ["lower_lim", "scale", "shift", "upper_lim"]
        self.tombo_read_start_rel_to_raw = "read_start_rel_to_raw"
        self.tombo_no = False;

    def m_initial(self):
        self.tombo_no = False;

    def read_signals(self, h5reader, scale_shift=None):
        self.m_initial();

        self.read_idx = list(h5reader[self.signal_full_path].keys())[0]
        raw_signals = h5reader[self.signal_full_path][self.read_idx][self.signal_entry][()]
        #print("Signals = {}".format(raw_signals[:10]))

        if scale_shift==None:
            channel_keyInfos = self.read_channel(h5reader)
            normal_signals = np.round((raw_signals+channel_keyInfos[self.channel_keyKeys[1]])/(channel_keyInfos[self.channel_keyKeys[0]]/channel_keyInfos[self.channel_keyKeys[2]]), 5)
        else:
            normal_signals = np.round((raw_signals-scale_shift[1])/scale_shift[0], 5)

        return (raw_signals, normal_signals)
    def read_fq(self, h5reader):
        self.m_initial();

        fq_4 = h5reader[self.fq_full_path][()].decode('UTF-8').splitlines()
        #print(fq_4)
        return fq_4
    def read_tombo_event(self, h5reader):
        self.m_initial();
        if self.tombo_event_full_path not in h5reader:
            self.tombo_no = True;
            return ;

        tombo_events = h5reader[self.tombo_event_full_path][:]
        '''print(type(tombo_events))
        print(tombo_events.shape)
        print(tombo_events)

        align_seq = []
        for _tei in range(len(tombo_events)):
            align_seq.append( tombo_events[_tei][4].decode('UTF-8') )
        print("Align_seq = {}".format(''.join(align_seq)))'''

        self.read_start_rel_to_raw = h5reader[self.tombo_event_full_path].attrs[self.tombo_read_start_rel_to_raw]

        for _tei in range(len(tombo_events)):
            tombo_events[_tei][2] += self.read_start_rel_to_raw

        return tombo_events
    def read_tombo_info(self, h5reader):
        self.m_initial();
        if self.tombo_align_full_path not in h5reader:
            self.tombo_no = True;
            return ;

        tombo_align_info = h5reader[self.tombo_align_full_path].attrs
        tombo_align_summary = []
        for _taa in self.tombo_align_attr:
            tombo_align_summary.append( tombo_align_info[_taa] )

        tombo_scale_info = h5reader[self.tombo_base_full_path].attrs
        tombo_scale_summary = []
        for _tss in self.tombo_signal_norm_keys:
            tombo_scale_summary.append( float(tombo_scale_info[_tss]) )
        #print(tombo_align_summary)
        return (tombo_align_summary, tombo_scale_summary)

    def read_channel(self, h5reader):
        self.m_initial();

        channel_keyInfos = h5reader[self.channel_path].attrs
        '''print(channel_keyInfos[self.channel_keyKeys[0]], channel_keyInfos[self.channel_keyKeys[1]],
          channel_keyInfos[self.channel_keyKeys[2]], channel_keyInfos[self.channel_keyKeys[3]],
          channel_keyInfos[self.channel_keyKeys[4]].decode("UTF-8"), type(channel_keyInfos[self.channel_keyKeys[4]]))'''
        return channel_keyInfos

class SeqReverse:
    def __init__(self):
        self.DNA_bp = {'A': 'T', 'C': 'G',
                       'T': 'A', 'G': 'C' }
        self.RNA_bp = {'A': 'U', 'C': 'G',
                       'U': 'A', 'G': 'C'}

        self.dna_bp_to_int_dict = {'A':0, 'C':1, 'G':2, 'T':3, \
                              'a':0, 'c':1, 'g':2, 't':3}
        self.rna_bp_to_int_dict = {'A':0, 'C':1, 'G':2, 'U':3, \
                              'a':0, 'c':1, 'g':2, 'u':3}

        self.rna_bp_to_int_dict_for_tombo = {'A':0, 'C':1, 'G':2, 'T':3, \
                                             'a':0, 'c':1, 'g':2, 't':3}

    def get_base_ind(self, mbase, mstrand, DNA_or_RNA, is_tombo=True):
        if DNA_or_RNA=='DNA':
            return self.dna_bp_to_int_dict[ mbase if mstrand=='+' else self.DNA_bp[mbase] ]
        else:
            if is_tombo:
                return self.rna_bp_to_int_dict_for_tombo[ mbase if mstrand=='+' else self.DNA_bp[mbase] ]
            else:
                return self.rna_bp_to_int_dict[ mbase if mstrand=='+' else self.RNA_bp[mbase] ]

    def get_dna_rev(self, seq):
        rseq = []
        for _i in range(len(seq)):
            rseq.append(self.DNA_bp[seq[_i]])
        return ''.join(rseq[::-1])
    def get_rna_rev(self, seq):
        rseq = []
        for _i in range(len(seq)):
            rseq.append(self.RNA_bp[seq[_i]])
        return ''.join(rseq[::-1])

'''
            DATASET "Events" {
               DATATYPE  H5T_COMPOUND {
                  H5T_IEEE_F64LE "norm_mean";
                  H5T_IEEE_F64LE "norm_stdev";
                  H5T_STD_U32LE "start";
                  H5T_STD_U32LE "length";
                  H5T_STRING {
                     STRSIZE 1;
                     STRPAD H5T_STR_NULLPAD;
                     CSET H5T_CSET_ASCII;
                     CTYPE H5T_C_S1;
                  } "base";
               }
               DATASPACE  SIMPLE { ( 248 ) / ( 248 ) }
               DATA {
               (0): {
                     -0.304124,
                     nan,
                     0,
                     51,
                     "A"
                  },

'''

if __name__ == '__main__':
    f5folder = sys.argv[1]
    f5r = Fast5Read()
    seqr = SeqReverse()

    sc_ref_file = '/data/data1/NanoporeRNA/rna_data.ref77/sc_rna_coding.fa'
    sc_ref = read_fasta.read_ref(sc_ref_file)

    noinlist =[];
    f5files = glob.glob(os.path.join(f5folder, "*/*.fast5"));
    for f5f in f5files:
        with h5py.File(f5f, 'r') as h5reader:
            fq_4 = f5r.read_fq( h5reader )
            tombo_sum = f5r.read_tombo_info( h5reader )
            if f5r.tombo_no: continue;
            tombo_align_summary =tombo_sum[0]
            tombo_scale_summary = tombo_sum[1]

            if tombo_align_summary[0]>5 or tombo_align_summary[1]>5: continue;
            if tombo_align_summary[5] not in noinlist:
                print(f5f, tombo_align_summary)
                tombo_events = f5r.read_tombo_event( h5reader)
                align_seq = []
                for _tei in range(len(tombo_events)):
                    align_seq.append( tombo_events[_tei][4].decode('UTF-8') )
                print("Align_seq = {}".format(''.join(align_seq)))
                if tombo_align_summary[5] =='+': print("Nanop_seq = {}".format(fq_4[1].replace('U','T')[tombo_align_summary[1]:(len(fq_4[1])-tombo_align_summary[0])]))
                else: print("Nanop_seq = {}".format(fq_4[1].replace('U','T')[tombo_align_summary[0]:(len(fq_4[1])-tombo_align_summary[1])]))
                #print("Nanop_seq = {}".format(fq_4[1].replace('U','T'))); #[tombo_align_summary[1]:(-tombo_align_summary[0])]))
                if tombo_align_summary[5] =='+': print("Refer_seq = {}".format(sc_ref[tombo_align_summary[2]][tombo_align_summary[4]:tombo_align_summary[3]]))
                else: print("Refer_seq = {}".format(seqr.get_dna_rev(sc_ref[tombo_align_summary[2]][tombo_align_summary[4]:tombo_align_summary[3]])))
                print((''.join(align_seq))==(sc_ref[tombo_align_summary[2]][tombo_align_summary[4]:tombo_align_summary[3]] if tombo_align_summary[5] =='+' else seqr.get_dna_rev(sc_ref[tombo_align_summary[2]][tombo_align_summary[4]:tombo_align_summary[3]]) ) )
                print(len(fq_4[1].replace('U','T')[tombo_align_summary[1]:(len(fq_4[1])-tombo_align_summary[0])]), len(align_seq))
            
                print(tombo_events[:10])
                noinlist.append(tombo_align_summary[5])
                if len(noinlist)==2: break;

