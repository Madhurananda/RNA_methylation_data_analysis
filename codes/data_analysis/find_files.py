
import os;
import sys;
sys.path.insert(0, '/home/madhu/work/codes')

from main import check_endWith, check_create_dir, check_env, write_file, count_lines
from tqdm import tqdm
import glob

def get_N_fast5(data_folder):
    all_fast5_files = [os.path.join(root, name)
        for root, dirs, files in os.walk(data_folder)
        for name in files
        if name.endswith((".fast5"))]

    print('Len of all_fast5_files: ', len(all_fast5_files), ' for: ', data_folder)

    return all_fast5_files




# fast5_files = glob.glob('/home/madhu/work/datasets/p3_m6A_RNA_modification_native_RNA_seq/*/extracted_data/*/.fast5')

# folder_to_find = ''



'''
This script will generate the final summary of the BAM files, Tombo .... 
'''

## Paper 3 

'/home/madhu/work/datasets/p3_m6A_RNA_modification_native_RNA_seq/GSM3528749/extracted_data/'

extracteddata_folder='/home/madhu/work/datasets/p3_m6A_RNA_modification_native_RNA_seq/'

# all_fast5_files = get_N_fast5(extracteddata_folder)



a_dataset = 'GSM3528749/extracted_data'
data_folder=extracteddata_folder+a_dataset
data_folder = check_endWith(data_folder)
all_fast5_files = get_N_fast5(data_folder)

# filesToFind = ['GXB01170_20180726_FAH85615_GA10000_sequencing_run_RNAAB089716_65817_read_14588_ch_350_strand', 'GXB01170_20180726_FAH85615_GA10000_sequencing_run_RNAAB089716_65817_read_15946_ch_375_strand', 'GXB01170_20180726_FAH85615_GA10000_sequencing_run_RNAAB089716_65817_read_20996_ch_257_strand', 'GXB01170_20180726_FAH85615_GA10000_sequencing_run_RNAAB089716_65817_read_21936_ch_283_strand', 'GXB01170_20180726_FAH85615_GA10000_sequencing_run_RNAAB089716_65817_read_22777_ch_458_strand', 'GXB01170_20180726_FAH85615_GA10000_sequencing_run_RNAAB089716_65817_read_17707_ch_277_strand', 'GXB01170_20180726_FAH85615_GA10000_sequencing_run_RNAAB089716_65817_read_26741_ch_313_strand', 'GXB01170_20180726_FAH85615_GA10000_sequencing_run_RNAAB089716_65817_read_18370_ch_281_strand', 'GXB01170_20180726_FAH85615_GA10000_sequencing_run_RNAAB089716_65817_read_17423_ch_143_strand', 'GXB01170_20180726_FAH85615_GA10000_sequencing_run_RNAAB089716_65817_read_20024_ch_253_strand']





a_dataset = 'GSM3528750/extracted_data'
data_folder=extracteddata_folder+a_dataset
data_folder = check_endWith(data_folder)
all_fast5_files.extend(get_N_fast5(data_folder))

filesToFind = ['GXB01170_20180726_FAH84534_GA20000_sequencing_run_RNAAB090763_92733_read_10152_ch_155_strand', 'GXB01170_20180726_FAH84534_GA20000_sequencing_run_RNAAB090763_92733_read_11002_ch_135_strand', 'GXB01170_20180726_FAH84534_GA20000_sequencing_run_RNAAB090763_92733_read_9599_ch_396_strand', 'GXB01170_20180726_FAH84534_GA20000_sequencing_run_RNAAB090763_92733_read_9768_ch_242_strand', 'GXB01170_20180726_FAH84534_GA20000_sequencing_run_RNAAB090763_92733_read_9188_ch_203_strand', 'GXB01170_20180726_FAH84534_GA20000_sequencing_run_RNAAB090763_92733_read_10830_ch_232_strand', 'GXB01170_20180726_FAH84534_GA20000_sequencing_run_RNAAB090763_92733_read_10543_ch_93_strand', 'GXB01170_20180726_FAH84534_GA20000_sequencing_run_RNAAB090763_92733_read_9067_ch_368_strand', 'GXB01170_20180726_FAH84534_GA20000_sequencing_run_RNAAB090763_92733_read_11467_ch_44_strand', 'GXB01170_20180726_FAH84534_GA20000_sequencing_run_RNAAB090763_92733_read_11521_ch_442_strand']






a_dataset = 'GSM3528751/extracted_data'
data_folder=extracteddata_folder+a_dataset
data_folder = check_endWith(data_folder)
all_fast5_files.extend(get_N_fast5(data_folder))

# filesToFind = ['MinION_2_20181108_FAK35406_MN29046_sequencing_run_RNA081120181_29618_read_30540_ch_494_strand.fast5', 'MinION_2_20181108_FAK35406_MN29046_sequencing_run_RNA081120181_29618_read_34597_ch_101_strand.fast5', 'MinION_2_20181108_FAK35406_MN29046_sequencing_run_RNA081120181_29618_read_38882_ch_379_strand', 'MinION_2_20181108_FAK35406_MN29046_sequencing_run_RNA081120181_29618_read_20916_ch_263_strand', 'MinION_2_20181108_FAK35406_MN29046_sequencing_run_RNA081120181_29618_read_38262_ch_369_strand', 'MinION_2_20181108_FAK35406_MN29046_sequencing_run_RNA081120181_29618_read_28602_ch_410_strand', 'MinION_2_20181108_FAK35406_MN29046_sequencing_run_RNA081120181_29618_read_33094_ch_157_strand', 'MinION_2_20181108_FAK35406_MN29046_sequencing_run_RNA081120181_29618_read_31399_ch_319_strand', 'MinION_2_20181108_FAK35406_MN29046_sequencing_run_RNA081120181_29618_read_21041_ch_263_strand', 'MinION_2_20181108_FAK35406_MN29046_sequencing_run_RNA081120181_29618_read_30294_ch_5_strand']







a_dataset = 'GSM3528752/extracted_data'
data_folder=extracteddata_folder+a_dataset
data_folder = check_endWith(data_folder)
all_fast5_files.extend(get_N_fast5(data_folder))

# filesToFind = ['MinION_2_20181108_FAK35429_MN29580_sequencing_run_RNA081120182_38847_read_39309_ch_12_strand', 'MinION_2_20181108_FAK35429_MN29580_sequencing_run_RNA081120182_38847_read_52717_ch_216_strand', 'MinION_2_20181108_FAK35429_MN29580_sequencing_run_RNA081120182_38847_read_32121_ch_461_strand', 'MinION_2_20181108_FAK35429_MN29580_sequencing_run_RNA081120182_38847_read_38389_ch_78_strand', 'MinION_2_20181108_FAK35429_MN29580_sequencing_run_RNA081120182_38847_read_54779_ch_432_strand', 'MinION_2_20181108_FAK35429_MN29580_sequencing_run_RNA081120182_38847_read_46917_ch_193_strand', 'MinION_2_20181108_FAK35429_MN29580_sequencing_run_RNA081120182_38847_read_46364_ch_344_strand', 'MinION_2_20181108_FAK35429_MN29580_sequencing_run_RNA081120182_38847_read_40516_ch_4_strand', 'MinION_2_20181108_FAK35429_MN29580_sequencing_run_RNA081120182_38847_read_36948_ch_5_strand', 'MinION_2_20181108_FAK35429_MN29580_sequencing_run_RNA081120182_38847_read_39707_ch_336_strand']



for a_file in filesToFind: 
    for an_exist_file in (all_fast5_files):
        if a_file in an_exist_file:
            print('an_exist_file: ', an_exist_file)

## Number of FAST5 and FASTQ files are at: 
'/home/madhu/work/basecall/guppy_v5.0.11/p3_m6A_RNA_modification_native_RNA_seq/GSM3528749/tombo_resquiggle_check_output.txt'


## The base mapping rate can be found at: 
'less /home/madhu/work/analysis/p3_m6A_RNA_modification_native_RNA_seq/GSM3528752/GSM3528752.bam.stat'


# print('The number of single read files is: ')



## Paper 5 


