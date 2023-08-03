

import os, sys
sys.path.insert(0, '/home/madhu/work/codes')
from main import check_endWith, check_create_dir, check_env, read_file, write_file



'''
Make BED format for reference sequence
'''

if __name__=='__main__':

    '''
    I will execute something like this:

    python /home/madhu/work/codes/data_analysis/make_BED.py /home/madhu/work/ref_transcriptome/OTHERS/synthetic_construct/Best_Practices_dRNAseq_analysis/reference_fasta/curlcake_constructs.fasta

    python /home/madhu/work/codes/data_analysis/make_BED.py /home/madhu/work/ref_transcriptome/IVT_seq/IVT_seq.fa

    '''


    ref_seq = sys.argv[1]
    # ref_seq = '/home/madhu/work/ref_transcriptome/IVT_seq/IVT_seq.fa'

    BED_type=''
    if len(sys.argv) > 2:
        BED_type = 'DRUMMER'

    file_cnt = read_file(ref_seq)


    BED_file = ''
    
    for a_chunk in file_cnt.split('>'):
        if not a_chunk == '':
            chr = a_chunk.split('\n')[0]
            ncl_seq = a_chunk.split('\n')[1]
            # print('\n*** Chromosome: ', chr)
            # print('\n*** Nucleotide sequence: ', ncl_seq)
            # print('\n*** Checking the final chars: ', ncl_seq[-1])
            if BED_type == 'DRUMMER':
                print('I need to update BED_file_output_name later if I need to ... ')
                BED_file += chr+'\n'
                # BED_file_output_name = ref_seq.replace(ref_seq.split('/')[-1], ref_seq.split('/')[-1].replace( ref_seq.split('/')[-1].split('.')[-1], 'BED'))
            elif BED_type == '':
                BED_file += chr+'\t'+'0'+'\t'+str(len(ncl_seq)-1)+'\t'+chr+'\t'+'.'+'\t'+'.'+'\n'
                BED_file_output_name = ref_seq.replace(ref_seq.split('/')[-1], ref_seq.split('/')[-1].replace( ref_seq.split('/')[-1].split('.')[-1], 'BED'))



    write_file(BED_file_output_name, BED_file)

    print('BED file has been created at: ', BED_file_output_name)
    print('CONTENT: ------- ')
    print(BED_file)

