
import os;
import sys;
sys.path.insert(0, '/home/madhu/work/codes')

from main import check_endWith, check_create_dir, check_env, write_file, count_lines


def get_N_fast5(data_folder):
    all_fast5_files = [os.path.join(root, name)
        for root, dirs, files in os.walk(data_folder)
        for name in files
        if name.endswith((".fast5"))]

    print('Len of all_fast5_files: ', len(all_fast5_files), ' for: ', a_dataset)





'''
This script will generate the final summary of the BAM files, Tombo .... 
'''

## Paper 3 

'/home/madhu/work/datasets/p3_m6A_RNA_modification_native_RNA_seq/GSM3528749/extracted_data/'

extracteddata_folder='/home/madhu/work/datasets/p3_m6A_RNA_modification_native_RNA_seq/'

a_dataset = 'GSM3528749'
data_folder=extracteddata_folder+a_dataset
data_folder = check_endWith(data_folder)
get_N_fast5(data_folder)


a_dataset = 'GSM3528750'
data_folder=extracteddata_folder+a_dataset
data_folder = check_endWith(data_folder)
get_N_fast5(data_folder)

a_dataset = 'GSM3528751'
data_folder=extracteddata_folder+a_dataset
data_folder = check_endWith(data_folder)
get_N_fast5(data_folder)

a_dataset = 'GSM3528752'
data_folder=extracteddata_folder+a_dataset
data_folder = check_endWith(data_folder)
get_N_fast5(data_folder)


## Number of FAST5 and FASTQ files are at: 
'/home/madhu/work/basecall/guppy_v5.0.11/p3_m6A_RNA_modification_native_RNA_seq/GSM3528749/tombo_resquiggle_check_output.txt'


## The base mapping rate can be found at: 
'less /home/madhu/work/analysis/p3_m6A_RNA_modification_native_RNA_seq/GSM3528752/GSM3528752.bam.stat'


print('The number of single read files is: ')



## Paper 5 


