
import os, sys
sys.path.insert(0, '/home/madhu/work/codes')

from main import check_endWith, check_create_dir, check_env, read_file

import pysam
from datetime import datetime
from tqdm import tqdm
from multiprocessing.pool import ThreadPool, Pool
from multiprocessing import cpu_count

import matplotlib.pyplot as plt
import numpy as np

'''
Here, I will be investigating ..... update ... 
'''



def read_bam(bamf_file, mq_thr = 20, length_thr = 1000):
    readid_list = {}

    n_read = 0
    n_good = 0
    n_bad = 0

    list_infer_len_thr = []
    list_query_algn_len = []
    list_MAPQ = []

    samfile = pysam.AlignmentFile(bamf_file, "rb")
    for read in samfile.fetch():
        n_read += 1
        
        infer_len_thr = read.infer_read_length()
        query_algn_len = read.query_alignment_length
        mapQ = read.mapping_quality

        if mapQ < mq_thr or infer_len_thr<length_thr or query_algn_len<length_thr:
            n_bad += 1
            pass
        else:
            if read.query_name not in readid_list:
                readid_list[ read.query_name ] = []
            # else:
            #     print('read.query_name', read.query_name)
            #     print( read.reference_name, read.reference_start, read.is_reverse, read.mapping_quality, read.query_alignment_length )
            
            readid_list[ read.query_name ].append( ( read.reference_name, read.reference_start, read.is_reverse, mapQ, query_algn_len) )
            n_good += 1

        #if len(readid_list)>10:
        #     break;
        list_infer_len_thr.append(infer_len_thr)
        list_query_algn_len.append(query_algn_len)
        list_MAPQ.append(mapQ)

    samfile.close()
    return readid_list, n_read, n_good, n_bad, list_infer_len_thr, list_query_algn_len, list_MAPQ


if __name__ == '__main__':

    check_env('ACONDA')

    startTime = datetime.now()
    current_time = startTime.strftime("%Y/%m/%d at %H:%M:%S")
    print('\n\nThe script starting at: ' + str(current_time), ' \n\n' )


    analysis_folder = sys.argv[1]
    analysis_folder = check_endWith(analysis_folder)
    
    given_dataset_list = sys.argv[2]
    if given_dataset_list == 'ALL_IN_DIR':
        print('All the folders inside the given folder will be considered.')
        dataset_list = next(os.walk(analysis_folder))[1]
    ## I can also make another option --exclude 
    else:
        dataset_list = given_dataset_list.split(',')
        dataset_list = list(filter(None, dataset_list)) # remove empty strings 
    
    ref_genome = sys.argv[3]
    ref_file = read_file(ref_genome)
    # print(ref_file.split('\n'))

    LEN_THR = int(sys.argv[4])
    MQ_THR = int(sys.argv[5])

    list_chr = []
    for a_line in ref_file.split('\n'):
        if a_line.startswith('>'):
            list_chr.append(a_line.replace('>', ''))
    

    for a_dataset in tqdm(dataset_list):
        data_analysis_folder=analysis_folder+a_dataset
        data_analysis_folder = check_endWith(data_analysis_folder)
        
        bam_file = data_analysis_folder+a_dataset+'.bam'

        bam_analysis_dir = data_analysis_folder+'BAM_analysis/'
        check_create_dir(bam_analysis_dir)

        name_fig_infer_len_thr = bam_analysis_dir+a_dataset+'__infer_len_thr.png'
        name_fig_query_algn_len = bam_analysis_dir+a_dataset+'__query_algn_len.png'
        # name_fig_both = bam_analysis_dir+a_dataset+'__both.png'
        name_fig_both = bam_analysis_dir+a_dataset+'__both_hist.png'

        print('\n\n******** The BAM file is: ', bam_file)
        # bamfile = sys.argv[1]
        # ref_genome = '/mnt/labshare/share/reference_genome/EpiNano_Reference_sequences/cc.fasta'


        # print(list_chr)


        # bam_file = '/home/madhu/work/analysis/p3_m6A_RNA_modification_native_RNA_seq/GSM3528749/GSM3528749.bam'

        testbam, n_read, n_good, n_bad, list_infer_len_thr, list_query_algn_len, list_MAPQ = read_bam(bam_file, length_thr=LEN_THR, mq_thr=MQ_THR)

        print('The number of reads: ', n_read)
        print('The number of GOOD reads: ', n_good)
        print('The number of BAD reads: ', n_bad)
        print('The GOOD percentage: ', round((n_good/n_read), 4)*100, '%')

        print('infer_len_thr: ', len(list_infer_len_thr))
        print('query_algn_len: ', len(list_query_algn_len))

        # # plotting labelled histogram
        # plt.hist(list_query_algn_len, bins=5, edgecolor='black')
        # plt.xlabel('Query Align Length')
        # plt.ylabel('Count')
        # plt.show()
        # # Save plot
        # plt.show()
        # plt.savefig(name_fig_query_algn_len)
        # plt.close()

        # # plotting labelled histogram
        # plt.plot(list_infer_len_thr, 'r')
        # plt.plot(list_query_algn_len, 'b')
        # # plt.xlabel('Query Align Length')
        # plt.ylabel('Lengths')
        # plt.show()
        # # Save plot
        # plt.show()
        # plt.savefig(name_fig_both)
        # plt.close()



        # plotting length and mapping quality in a seperate subplots ... 

        # plt.plot(list_infer_len_thr, 'r')
        # plt.plot(list_query_algn_len, 'b')
        # # plt.xlabel('Query Align Length')
        # plt.ylabel('Lengths')
        # plt.show()
        # # Save plot
        # plt.show()
        # plt.savefig(name_fig_both)
        # plt.close()

        # fig, ax = plt.subplots(figsize=(20, 10))
        # plt.subplot(2, 1, 1)
        # plt.plot(list_query_algn_len)
        # plt.title('Mapping Quality and Lengths')
        # plt.ylabel('query_alignment_length')

        # plt.subplot(2, 1, 2)
        # plt.plot(list_MAPQ)
        # plt.xlabel('No. of Reads')
        # plt.ylabel('Maping quality')
        # plt.show()

        # plt.savefig(name_fig_both)
        # plt.close()

        fig, ax = plt.subplots(figsize=(20, 10))
        plt.subplot(1, 2, 1)

        # # fine ... 
        # plt.hist(list_query_algn_len, bins=2, range=[0, 2*LEN_THR], edgecolor='black')

        # list_BIN = [0, LEN_THR] + list(range(LEN_THR, max(list_query_algn_len), int((max(list_query_algn_len)-LEN_THR)/5))) +  [max(list_query_algn_len)]
        max_c_X_lim = 2000
        if max(list_query_algn_len) > max_c_X_lim:
            list_BIN_range = [0, LEN_THR] + [500, 1000, 1500, max_c_X_lim] + [max(list_query_algn_len)]
        else:
            list_BIN_range = [0, LEN_THR] + [500, 1000, 1500, max_c_X_lim] 
        print('list_BIN_range: ', list_BIN_range)
        counts, edges, bars = plt.hist(list_query_algn_len, bins=list_BIN_range, edgecolor='black')
        plt.bar_label(bars, fontsize=12)

        # print('counts: ', max(counts))

        # list(range(LEN_THR, max(list_query_algn_len), int((max(list_query_algn_len)-LEN_THR)/5)))

        # plt.hist(list_query_algn_len, bins=[0, LEN_THR, max(list_query_algn_len)], edgecolor='black')
        # plt.xlim(xmin=0, xmax = 2*LEN_THR)

        plt.xticks(list_BIN_range, list_BIN_range, fontsize=15)
        plt.yticks(np.arange(0,int(max(counts)),int(max(counts)/10)), np.arange(0,int(max(counts)),int(max(counts)/10)), fontsize=15)
        # plt.xticks(np.arange(len(list_BIN)), list_BIN)
        # plt.yticks(np.arange(0, 81, 10))
        # plt.legend((p1[0], p2[0]), ('Men', 'Women'))

        plt.xlabel('Query Align Length', fontsize=20)
        plt.ylabel('Counts', fontsize=20)

        x = [LEN_THR, LEN_THR]
        y = [0, max(counts)]
        plt.plot(x, y, color='red', linestyle='--')

        # plt.plot([LEN_THR, 0], [LEN_THR, 40000], color='red', linestyle='--')
        # plt.plot(list_query_algn_len)
        # plt.title('Mapping Quality and Lengths')
        # plt.ylabel('query_alignment_length')





        plt.subplot(1, 2, 2)


        max_c_X_lim = 80
        if max(list_MAPQ) > max_c_X_lim:
            list_BIN_range = [0, MQ_THR] + [40, 60, max_c_X_lim] + [max(list_MAPQ)]
        else:
            list_BIN_range = [0, MQ_THR] + [40, 60, max_c_X_lim] 
        print('list_BIN_range: ', list_BIN_range)
        counts, edges, bars = plt.hist(list_MAPQ, bins=list_BIN_range, edgecolor='black')
        plt.bar_label(bars, fontsize=12)

        plt.xticks(list_BIN_range, list_BIN_range, fontsize=15)
        plt.yticks(np.arange(0,int(max(counts)),int(max(counts)/10)), np.arange(0,int(max(counts)),int(max(counts)/10)), fontsize=15)
        # plt.xticks(np.arange(len(list_BIN)), list_BIN)
        # plt.yticks(np.arange(0, 81, 10))
        # plt.legend((p1[0], p2[0]), ('Men', 'Women'))

        plt.xlabel('Maping quality', fontsize=20)
        plt.ylabel('Counts', fontsize=20)

        x = [MQ_THR, MQ_THR]
        y = [0, max(counts)]
        plt.plot(x, y, color='red', linestyle='--')



        # plt.hist(list_MAPQ, bins=[0, MQ_THR, max(list_MAPQ)], edgecolor='black')
        # plt.xlabel('Maping quality', fontsize=20)
        # plt.ylabel('Counts', fontsize=20)
        # plt.plot(list_MAPQ)
        # plt.xlabel('No. of Reads')
        # plt.ylabel('Maping quality')
        plt.show()

        plt.savefig(name_fig_both)
        plt.close()



        chr_dict = {}
        for a_chr in list_chr:
            chr_dict[a_chr] = []
            for _rk in testbam.keys():
                if a_chr == testbam[_rk][0][0]:
                    chr_dict[a_chr].append(testbam[_rk][0][0])
                    # print(_rk, testbam[_rk]);
            print( 'The chr: ', a_chr, ' contains', len(chr_dict[a_chr]), ' reads.' )
        
        print('******** The end of the BAM file: ', bam_file)

    executionTime = (datetime.now() - startTime)
    current_time = datetime.now().strftime("%Y/%m/%d at %H:%M:%S")
    print('\n\nThe script completed at: ' + str(current_time))
    print('Execution time: ' + str(executionTime), ' \n\n')

