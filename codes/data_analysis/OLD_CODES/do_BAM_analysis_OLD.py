
import os
import sys
sys.path.insert(0, '/home/madhu/datasets/codes')

from main import check_endWith, check_create_dir, check_conda_env
from datetime import datetime
from tqdm import tqdm

'''
This script creates the BAM file, analyses it and produces statistics ... An example is given below which I had to do manually ...  


conda activate ACONDA

1. 
cat /home/madhu/datasets/basecall/guppy_v5.0.11/p2_identification_m6A_m5C_RNA_modification/HeLa_all_fast5_M11B_KO/*/*.fastq > /home/madhu/datasets/analysis/p2_identification_m6A_m5C_RNA_modification/HeLa_all_fast5_M11B_KO/HeLa_all_fast5_M11B_KO.fq


2. 
minimap2 -t 5 --cs -ax splice -uf -k14 /home/madhu/datasets/ref_transcriptome/human/hg38.fa /home/madhu/datasets/analysis/p2_identification_m6A_m5C_RNA_modification/HeLa_all_fast5_M11B_KO/HeLa_all_fast5_M11B_KO.fq | samtools sort -t 5 -o /home/madhu/datasets/analysis/p2_identification_m6A_m5C_RNA_modification/HeLa_all_fast5_M11B_KO/HeLa_all_fast5_M11B_KO.bam


3. 
samtools index /home/madhu/datasets/analysis/p2_identification_m6A_m5C_RNA_modification/HeLa_all_fast5_M11B_KO/HeLa_all_fast5_M11B_KO.bam


4. 
samtools view /home/madhu/datasets/analysis/p2_identification_m6A_m5C_RNA_modification/HeLa_all_fast5_M11B_KO/HeLa_all_fast5_M11B_KO.bam | less


5. 
samtools stats /home/madhu/datasets/analysis/p2_identification_m6A_m5C_RNA_modification/HeLa_all_fast5_M11B_KO/HeLa_all_fast5_M11B_KO.bam > /home/madhu/datasets/analysis/p2_identification_m6A_m5C_RNA_modification/HeLa_all_fast5_M11B_KO/HeLa_all_fast5_M11B_KO.bam.stat


6.
less /home/madhu/datasets/analysis/p2_identification_m6A_m5C_RNA_modification/HeLa_all_fast5_M11B_KO/HeLa_all_fast5_M11B_KO.bam.stat


7. conda activate lrst_py39
python /home/madhu/datasets/scripts/LongReadSum bam -i /home/madhu/datasets/analysis/p2_identification_m6A_m5C_RNA_modification/HeLa_all_fast5_M11B_KO/HeLa_all_fast5_M11B_KO.bam -o /home/madhu/datasets/analysis/p2_identification_m6A_m5C_RNA_modification/HeLa_all_fast5_M11B_KO/long_output/



'''


if __name__=='__main__':
    if len(sys.argv)<2:
        print("Usage: {} {}".format(sys.argv[0], "basefolder"))
        print("Example:")
        print("       python {} {} {} {}".format(sys.argv[0], "/home/madhu/datasets/download/m6A_RNA_modification_native_RNA_seq/GSM3528751/extracted_data"))
#           print("       python {} {} {} {}".format(sys.argv[0], "Curlcake_non"))
        sys.exit(1)
    
    check_conda_env('ACONDA')

    startTime = datetime.now()
    current_time = startTime.strftime("%Y/%m/%d at %H:%M:%S")
    print('\n\nThe script starting at: ' + str(current_time), ' \n\n' )

    ## Basecall folder should be the folder which has been used as the output of guppy. 
    # basecalled_folder = '/home/madhu/datasets/basecall/guppy_v5.0.11/p2_identification_m6A_m5C_RNA_modification/HeLa_all_fast5_M11A_KO'
    basecalled_folder = sys.argv[1]
    basecalled_folder = check_endWith(basecalled_folder)
    
    given_dataset_list = sys.argv[2]
    if given_dataset_list == 'ALL_IN_DIR':
        print('All the folders inside the given folder will be considered.')
        dataset_list = next(os.walk(basecalled_folder))[1]
    ## I can also make another option --exclude 
    else:
        dataset_list = given_dataset_list.split(',')
        dataset_list = list(filter(None, dataset_list)) # remove empty strings 
        # print(dataset_list)

    analysis_folder = sys.argv[3]
    analysis_folder = check_endWith(analysis_folder)

    ref_genome = sys.argv[4]

    # print(list(dataset_list))

    if not type(dataset_list) == list:
        print('Please input a list for the datasets.')
        print('Example: python do_BAM_analysis.py /home/madhu/datasets/basecall/guppy_v5.0.11/p2_identification_m6A_m5C_RNA_modification [HeLa_all_fast5_M11A_KO, HeLa_all_fast5_M11C_WT, HeLa_all_fast5_M11A_WT, HeLa_all_fast5_M11C_KO] /home/madhu/datasets/analysis/p2_identification_m6A_m5C_RNA_modification /home/madhu/datasets/ref_transcriptome/human/hg38.fa')
        print('Exiting .... ')
        sys.exit(1)

    for a_dataset in tqdm(dataset_list):
        
        ## analysis folder is the folder where all the BAM and other files will be ... 
        # analysis_folder = '/home/madhu/datasets/analysis/p2_identification_m6A_m5C_RNA_modification/HeLa_all_fast5_M11A_KO'
        
        # analysis_folder = check_endWith(analysis_folder)

        # Make the same folder in the analysis directories  
        data_analysis_folder=analysis_folder+a_dataset
        data_analysis_folder = check_endWith(data_analysis_folder)
        
        # print('The analysis folder is: ', data_analysis_folder)

        # create the folder if necessary
        check_create_dir(data_analysis_folder)


        fq_file = data_analysis_folder+a_dataset+'.fq'
        bam_file = data_analysis_folder+a_dataset+'.bam'
        sam_index_file = bam_file+'.bai'
        sam_stat_file = bam_file+'.stat'
        # print('fq_file: ', fq_file)


        data_basecalled_folder= basecalled_folder + a_dataset
        data_basecalled_folder = check_endWith(data_basecalled_folder)


        ## Step 1: create the list of fastq files. 
        if os.path.isfile(fq_file):
            print('\n\n'+fq_file + ' ALREADY EXISTS. ')
        else:
            cat_fastq_cmd = "cat {} > {}".format( data_basecalled_folder+'*/*.fastq', fq_file )
            os.system(cat_fastq_cmd)
            print('cat_fastq_cmd is: ', cat_fastq_cmd)
            print('\n\n'+fq_file +' has been successfully created.')


        ## Step 2: Minimap tool
        if os.path.isfile(bam_file):
            print('\n\n'+bam_file + ' bam_file ALREADY EXISTS. ')
        else:
            minimap_cmd = 'minimap2 -t 5 --cs -ax splice -uf -k14 {} {} | samtools sort -t 5 -o {}'.format(ref_genome, fq_file, bam_file)
            os.system(minimap_cmd)
            print('minimap_cmd is: ', minimap_cmd)
            print('\n\n'+bam_file + ' has been successfully created: ')


        ## Step 3: SamTools index
        if os.path.isfile(sam_index_file):
            print('\n\n'+sam_index_file + ' ALREADY EXISTS. ')
        else:
            samtools_index_cmd = 'samtools index {}'.format(bam_file)
            os.system(samtools_index_cmd)
            print('samtools_index_cmd is: ', samtools_index_cmd)
            print('\n\n'+sam_index_file + ' has been successfully created: ')


        ## Step 4: SamTools Stats
        if os.path.isfile(sam_stat_file):
            print('\n\n'+sam_stat_file + ' ALREADY EXISTS. ')
        else:
            samtools_stat_cmd = 'samtools stats {} > {}'.format(bam_file, sam_stat_file)
            os.system(samtools_stat_cmd)
            print('samtools_stat_cmd is: ', samtools_stat_cmd)
            print('\n\n'+sam_stat_file + ' has been successfully created: ')


        ## Step 5: Read the BAM stat file and produce the analysis ... 
        view_BAM_file = os.popen('less {}'.format(sam_stat_file)).read()

        # print('\n\n\n BAM file: '+ view_BAM_file+ '\n\n\n')



        total_len = int(view_BAM_file.split('total length:')[1].split('bases mapped:')[0].split('\t')[1])
        bases_mapped = int(view_BAM_file.split('bases mapped:')[1].split('bases mapped (cigar):')[0].split('\t')[1])


        print('total_len: ', total_len)
        print('bases_mapped: ', bases_mapped)

        # print(view_BAM_file.split('bases mapped:')[1].split('bases mapped (cigar):')[0].split('\t')[1])

        print('\n\n\n****** The percentage of bases mapped, for '+data_analysis_folder+' is: ', round(((bases_mapped/total_len)*100),2), '%\n\n\n')


        ## Step 7: I am not doing this as I can't install LongReadSum to ACONDA environment. 
        # Use the new script: do_LongReadSum.py 

    executionTime = (datetime.now() - startTime)
    current_time = datetime.now().strftime("%Y/%m/%d at %H:%M:%S")
    print('\n\nThe script completed at: ' + str(current_time))
    print('Execution time: ' + str(executionTime), ' \n\n')
