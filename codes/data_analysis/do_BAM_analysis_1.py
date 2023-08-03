
import os
import sys
from main import check_endWith, check_create_dir

'''
This script creates the BAM file, analyses it and produces statistics ... An example is given below which I had to do manually ...  


conda activate lrst_py39

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

# def check_endWith(basecalled_folder):
#     if not basecalled_folder.endswith('/'):
#         basecalled_folder+='/'
#     return basecalled_folder

# def check_create_dir(folder):
#     if not os.path.isdir(folder):
#         os.makedirs(folder, exist_ok=True)
#         print('**\nThe Folder has been created: ', folder, '**\n')
#     else:
#         print('The folder already exists: ', folder)


if __name__=='__main__':
    if len(sys.argv)<2:
        print("Usage: {} {}".format(sys.argv[0], "basefolder"))
        print("Example:")
        print("       python {} {} {} {}".format(sys.argv[0], "/home/madhu/datasets/download/m6A_RNA_modification_native_RNA_seq/GSM3528751/extracted_data"))
#           print("       python {} {} {} {}".format(sys.argv[0], "Curlcake_non"))
        sys.exit(1)
    

    ## Basecall folder should be the folder which has been used as the output of guppy. 
    # basecalled_folder = '/home/madhu/datasets/basecall/guppy_v5.0.11/p2_identification_m6A_m5C_RNA_modification/HeLa_all_fast5_M11A_KO'
    basecalled_folder = sys.argv[1]
    basecalled_folder = check_endWith(basecalled_folder)


    ## analysis folder is the folder where all the BAM and other files will be ... 
    # analysis_folder = '/home/madhu/datasets/analysis/p2_identification_m6A_m5C_RNA_modification/HeLa_all_fast5_M11A_KO'
    analysis_folder = sys.argv[2]
    analysis_folder = check_endWith(analysis_folder)

    # Make the same folder in the analysis directories  
    analysis_folder+=basecalled_folder.split('/')[ len(basecalled_folder.split('/'))-2]
    analysis_folder = check_endWith(analysis_folder)
    
    # create the folder if necessary
    check_create_dir(analysis_folder)
    # if not os.path.isdir(analysis_folder):
    #     os.makedirs(analysis_folder, exist_ok=True)
    #     print('**\nThe Folder has been created: ', analysis_folder, '**\n')
    # else:
    #     print('The folder already exists: ', analysis_folder)


    # ref_genome = '/home/madhu/datasets/ref_transcriptome/human/hg38.fa'
    ref_genome = sys.argv[3]


    fq_file = analysis_folder+analysis_folder.split('/')[ len(analysis_folder.split('/'))-2]+'.fq'
    bam_file = analysis_folder+analysis_folder.split('/')[ len(analysis_folder.split('/'))-2]+'.bam'
    sam_index_file = bam_file+'.bai'
    sam_stat_file = bam_file+'.stat'
    # print('bam_file: ', bam_file)



    ## Step 1: create the list of fastq files. 

    # print(basecalled_folder+'*/*.fastq')
    # print(analysis_folder+analysis_folder.split('/')[ len(analysis_folder.split('/'))-2]+'.fq')

    if os.path.isfile(fq_file):
        print('\n\n'+fq_file + ' ALREADY EXISTS. ')
    else:
        os.system("cat {} > {}".format( basecalled_folder+'*/*.fastq', fq_file ));
        print('\n\n'+fq_file +' has been successfully created.')



    ## Step 2: Minimap tool
    if os.path.isfile(bam_file):
        print('\n\n'+bam_file + ' bam_file ALREADY EXISTS. ')
    else:
        minimap_cmd = 'minimap2 -t 5 --cs -ax splice -uf -k14 {} {} | samtools sort -t 5 -o {}'.format(ref_genome, fq_file, bam_file)
        os.system(minimap_cmd)
        print('\n\n'+bam_file + ' has been successfully created: ')


    ## Step 3: SamTools index
    if os.path.isfile(sam_index_file):
        print('\n\n'+sam_index_file + ' ALREADY EXISTS. ')
    else:
        samtools_index_cmd = 'samtools index {}'.format(bam_file)
        os.system(samtools_index_cmd)
        print('\n\n'+sam_index_file + ' has been successfully created: ')


    ## Step 4: SamTools Stats
    if os.path.isfile(sam_stat_file):
        print('\n\n'+sam_stat_file + ' ALREADY EXISTS. ')
    else:
        samtools_stat_cmd = 'samtools stats {} > {}'.format(bam_file, sam_stat_file)
        os.system(samtools_stat_cmd)
        print('\n\n'+sam_stat_file + ' has been successfully created: ')


    ## Step 5: Read the BAM stat file and produce the analysis ... 
    view_BAM_file = os.popen('less {}'.format(sam_stat_file)).read()

    # print('\n\n\n BAM file: '+ view_BAM_file+ '\n\n\n')



    total_len = int(view_BAM_file.split('total length:')[1].split('bases mapped:')[0].split('\t')[1])
    bases_mapped = int(view_BAM_file.split('bases mapped:')[1].split('bases mapped (cigar):')[0].split('\t')[1])


    print('total_len: ', total_len)
    print('bases_mapped: ', bases_mapped)

    # print(view_BAM_file.split('bases mapped:')[1].split('bases mapped (cigar):')[0].split('\t')[1])

    print('\n\n\n****** The percentage of bases mapped is: ', round(((bases_mapped/total_len)*100),4), '%\n\n\n')


    ## Step 7: I am not doing this as I can't install LongReadSum to ACONDA environment. 
    # Use the new script: do_LongReadSum.py 

