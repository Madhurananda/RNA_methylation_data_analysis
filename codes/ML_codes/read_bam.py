
import os, sys
import pysam

def read_bam(bamf_file, mq_thr = 20, length_thr = 1000):
    readid_list = {}

    samfile = pysam.AlignmentFile(bamf_file, "rb")
    for read in samfile.fetch():
        if read.mapping_quality < mq_thr or read.infer_read_length()<length_thr or read.query_alignment_length<length_thr:
            pass
        else:
            if read.query_name not in readid_list:
                readid_list[ read.query_name ] = []
            
            readid_list[ read.query_name ].append( ( read.reference_name, read.reference_start, read.is_reverse, read.mapping_quality, read.query_alignment_length) )

        #if len(readid_list)>10:
        #     break;
    samfile.close()
    return readid_list


if __name__ == '__main__':
    bamfile = sys.argv[1]

    testbam = read_bam(bamfile)

    for _rk in testbam.keys():
        print(_rk, testbam[_rk]);



