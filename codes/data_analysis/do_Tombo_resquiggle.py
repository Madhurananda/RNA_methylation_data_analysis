import os
import sys
sys.path.insert(0, '/home/madhu/work/codes')

from main import check_endWith, check_create_dir, check_env, write_file, count_lines
from datetime import datetime
from tqdm import tqdm
import h5py
from multiprocessing.pool import ThreadPool, Pool
from multiprocessing import cpu_count


'''
This script runs the Tombo resquiggle for each dataset for a given directory. 

conda activate ont4

(ont4) madhu@qgenlabgpu:~$ tombo resquiggle /home/madhu/datasets/basecall/guppy_v5.0.11/m6A_RNA_modification_native_RNA_seq/GSM3528752/ /mnt/labshare/share/reference_genome/EpiNano_Reference_sequences/cc.fasta --processes 10 --num-most-common-errors 5

inputs: 'python' 'do_Tombo_resquiggle.py' 'basecall folder where the FAST5 files are' 'all the dataset seperated by comma, no space' 'reference genome'

python /home/madhu/datasets/codes/data_analysis/do_Tombo_resquiggle.py /home/madhu/datasets/basecall/guppy_v5.0.11/p2_identification_m6A_m5C_RNA_modification HeLa_all_fast5_M11B_KO,HeLa_all_fast5_M11A_KO,HeLa_all_fast5_M11C_WT,HeLa_all_fast5_M11A_WT,HeLa_all_fast5_M11C_KO,HeLa_all_fast5_M11B_WT /home/madhu/datasets/ref_transcriptome/human/hg38.fa 20

tombo resquiggle /home/madhu/datasets/basecall/guppy_v5.0.11/p2_identification_m6A_m5C_RNA_modification/HeLa_all_fast5_M11B_KO/ /home/madhu/datasets/ref_transcriptome/human/hg38.fa --processes 10 --num-most-common-errors 5

'''



def do_Tombo(args):
    t0 = datetime.now()
    a_dataset, basecalled_folder, ref_genome, n_CPU = args[0], args[1], args[2], args[3]
    try:
        # Make the same folder in the analysis directories  
        data_folder=basecalled_folder+a_dataset
        data_folder = check_endWith(data_folder)
        
        to_wrtie_content = ''
        write_file_name = data_folder+'tombo_resquiggle_check_output.txt'
        
        Tombo_cmd = 'tombo resquiggle {} {} --processes {} --overwrite --num-most-common-errors 5'.format(data_folder, ref_genome, n_CPU)
        print('\nTombo_cmd :', Tombo_cmd)
        # os.system(Tombo_cmd)
        # print('\n\n'+Tombo_cmd + ' has been successfully executed for: ', data_folder)
        stream = os.popen(Tombo_cmd)
        output = stream.read()
        print('\n\nThe TOMBO output: ------------------------------------------------ \n', (output), '\n\n')
        to_wrtie_content += '\n\nThe TOMBO output: ------------------------------------------------ \n' + str(output) + '\n\n'
        sys.stdout.flush()



        '''
        Now, here also do the Tombo statistics ..., same as check_tombo_resquiggle.py .... 
        26th April, 2023: I decided to skip checking each FAST5 files as it is taking a long time and miltiprocessing inside another multitasking didn't work better. 
        '''
        
        '''
        Just comment out to make these steps work ... 
        '''
        # ## Step 1: Read all the fast5 files in the basefolder 
        # all_fast5_files = [os.path.join(root, name)
        #      for root, dirs, files in os.walk(data_folder)
        #      for name in files
        #      if name.endswith((".fast5"))]
        # print('Len of all_fast5_files: ', len(all_fast5_files))
        # to_wrtie_content += ('\nLen of all_fast5_files: '+str(len(all_fast5_files)) )
        # sys.stdout.flush()
        # ## Step 2: Calculate Tombo performance by analysing the FAST5 files 
        # n_tombo = 0
        # n_Ntombo =0
        # for a_fast5_file in tqdm(all_fast5_files):
        #     f5r = h5py.File(a_fast5_file, 'r')
        #     if '/Analyses/RawGenomeCorrected_000/BaseCalled_template/Events' in f5r:
        # #          print('This file is done with Tombo: ', a_fast5_file)
        #         n_tombo += 1
        #     else:
        # #          print('This file is NOT done with Tombo: ', a_fast5_file)
        #         n_Ntombo += 1
        
        # print('\n\n', a_dataset, ': The no of files with Tombo: ')
        # print('\n\n', a_dataset, ': The no of files with NO Tombo: ')
        
        # print('\n\n\n', a_dataset, ': The percentage of files for which tombo done: ', round(((n_tombo/(n_tombo+n_Ntombo))*100), 2), '%\n\n')
        
        # to_wrtie_content += ('\n'+a_dataset+': The no of files with Tombo: '+str(n_tombo) )
        # to_wrtie_content += ('\n'+a_dataset+': The no of files with NO Tombo: '+str(n_Ntombo) )
        # to_wrtie_content += ('\n\n'+a_dataset+': The percentage of files for which tombo done:  '+str(round(((n_tombo/(n_tombo+n_Ntombo))*100), 2)) +'%\n\n')
        # sys.stdout.flush()



        # # Step 3:  Get all the FASTQ files
        # all_fastQ_files = [os.path.join(root, name)
        #             for root, dirs, files in os.walk(data_folder)
        #             for name in files
        #             if name.endswith((".fastq"))]
        # print('Len of all_fastQ_files: ', len(all_fastQ_files))
        # to_wrtie_content += ('\nLen of all_fastQ_files: '+str(len(all_fastQ_files)) )
        # sys.stdout.flush()

        # # Step 4:  Get the number of lines from each of those FASTQ files
        # c_pass = 0
        # c_fail = 0
        # for a_fastQ_file in tqdm(all_fastQ_files): 
        # #       print('The fastQ file is: ', a_fastQ_file)
        #     if data_folder+'pass/' in a_fastQ_file:
        # #          print('The PASS fastQ file is: ', a_fastQ_file)
        #         n_lines = int(count_lines(a_fastQ_file))
        #         c_pass += n_lines
        #     elif data_folder+'fail/' in a_fastQ_file:
        # #          print('The FAIL fastQ file is: ', a_fastQ_file)
        #         n_lines = int(count_lines(a_fastQ_file))
        #         c_fail += n_lines
        #     else:
        #         print('\n\n***** There is something VERY wrong::', a_fastQ_file, ':: INVESTIGATE...******\n\n')
        # sys.stdout.flush()


        # n_read_pass = c_pass/4
        # n_read_fail = c_fail/4
        
        # print('\n', a_dataset, ': The total pass line counts are: ', c_pass)
        # print('\n', a_dataset, ': The total fail line counts are: ', c_fail)
        # print('\n\n', a_dataset, ': pass/(pass+fail) READ percentage is: ', round((n_read_pass/(n_read_pass+n_read_fail))*100, 2), '%')

        # to_wrtie_content += ('\n'+a_dataset+': The total pass line counts are: '+str(c_pass) )
        # to_wrtie_content += ('\n'+a_dataset+': The total fail line counts are: '+str(c_fail) )
        # to_wrtie_content += ('\n\n'+a_dataset+': pass/(pass+fail) READ percentage is: '+str(round((n_read_pass/(n_read_pass+n_read_fail))*100, 2)) + '%\n\n' )
        # write_file(write_file_name, to_wrtie_content)
        # sys.stdout.flush()

    except Exception as e:
        print('a_dataset: ', a_dataset)
        print('Exception:', e)
    
    sys.stdout.flush()
    return(a_dataset, datetime.now() - t0)





if __name__=='__main__':
    if len(sys.argv)<2:
        print("Usage: {} {}".format(sys.argv[0], "basefolder"))
        print("Example:")
        print("       python {} {} {} {}".format(sys.argv[0], "/home/madhu/datasets/download/m6A_RNA_modification_native_RNA_seq/GSM3528751/extracted_data"))
#           print("       python {} {} {} {}".format(sys.argv[0], "Curlcake_non"))
        sys.exit(1)
    
    check_env('ont4')

    startTime = datetime.now()
    current_time = startTime.strftime("%Y/%m/%d at %H:%M:%S")
    print('\n\nThe script starting at: ' + str(current_time), ' \n\n' )

    ## Basecall folder should be the folder which has been used as the output of guppy. 
    # basecalled_folder = '/home/madhu/datasets/basecall/guppy_v5.0.11/p2_identification_m6A_m5C_RNA_modification/HeLa_all_fast5_M11A_KO'
    basecalled_folder = sys.argv[1]
    basecalled_folder = check_endWith(basecalled_folder)

    # dataset_list = sys.argv[2].split(',')
    # dataset_list = list(filter(None, dataset_list)) # remove empty strings 
    # # print(dataset_list)

    given_dataset_list = sys.argv[2]
    if given_dataset_list == 'ALL_IN_DIR':
        print('All the folders inside the given folder will be considered.')
        dataset_list = next(os.walk(basecalled_folder))[1]
    ## I can also make another option --exclude 
    else:
        dataset_list = given_dataset_list.split(',')
        dataset_list = list(filter(None, dataset_list)) # remove empty strings 
        # print(dataset_list)

    ref_genome = sys.argv[3]
    n_CPU = sys.argv[4]

    inputs = zip(dataset_list, [basecalled_folder]*len(dataset_list), [ref_genome]*len(dataset_list), [n_CPU]*len(dataset_list))

    ## Use the maximum CPU possible 
    cpus = cpu_count()
    n_jobs = min(cpus, len(dataset_list))

    # n_jobs = 20

    print('n_jobs: ', n_jobs)
    ## I am using multiple Processors as it is faster than threads. 
    results = tqdm(Pool(n_jobs).imap_unordered(do_Tombo, inputs), total=len(dataset_list))
    # results = ThreadPool(n_jobs).imap_unordered(unzip, inputs)
    for result in results:
        # print('time (s):', result)
        print('File: ', result[0], 'time :', result[1])
    

    executionTime = (datetime.now() - startTime)
    current_time = datetime.now().strftime("%Y/%m/%d at %H:%M:%S")
    print('\n\nThe script completed at: ' + str(current_time))
    print('Execution time: ' + str(executionTime), ' \n\n')


