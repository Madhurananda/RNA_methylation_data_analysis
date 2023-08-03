
import os
import sys

sys.path.insert(0, '/home/madhu/work/codes')
from main import check_endWith, check_create_dir, check_env, write_file

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
from tqdm import tqdm
from operator import itemgetter
from datetime import datetime
from multiprocessing.pool import ThreadPool, Pool
from multiprocessing import cpu_count

'''
I am creating this file to cehck the MAPQ (mapping quality) in a BAM file. 
This file should produce a distribution of all the MAPQ values. 
'''

def do_MAPQ(args):
    t0 = datetime.now()
    BAM_file = args[0]
    try:
        output_dir = BAM_file.replace(BAM_file.split('/')[ len(BAM_file.split('/'))-1], '')
        name_fig = BAM_file.split('.bam')[0]+'_MAPQ.png'
        name_txt = BAM_file.split('.bam')[0]+'_MAPQ_info.txt'

        if os.path.exists(name_fig) and os.path.exists(name_txt):
            print('\n\n**** Checking has already been done for: ', output_dir)
        else:
            # output_dir = BAM_file.replace(BAM_file.split('/')[ len(BAM_file.split('/'))-1], 'long_output/') 
            print('BAM_file: ', BAM_file)
            print('out_dir: ', output_dir)
            print('\n')
            
            samtool_cmd = 'samtools view {} | cut -f5 | sort -n | uniq -c'.format(BAM_file)
            stream = os.popen(samtool_cmd)
            output = stream.read()
            print('samtool_cmd is: ', samtool_cmd)
            print('\n\nThe MAPQ distribution for the file: ', BAM_file, ' is: \n', (output))

            # Write the output to the file ... 
            to_wrtie_content = '\n\nThe MAPQ distribution for the file: '+ str(BAM_file)+ ' is: \n'+ str(output) + '\n\n'
            write_file(name_txt, to_wrtie_content)
            
            
            list_count = []
            list_mapQ = []

            for a_line in output.split('\n'):
                if a_line == '':
                    continue
                s_line = a_line.split(' ')
                while('' in s_line):
                    s_line.remove('')
            #     print('line: ', s_line)
                count = s_line[0]
                mapQ = s_line[1]
            #     print('count: ', count)
            #     print('mapQ: ', mapQ)
            #     print('\n\n')
                list_count.append(int(count))
                list_mapQ.append(int(mapQ))

            ## Now, as I have sorted out the array in the terminal, no need to sort out in the Python ... 
            x = list_mapQ
            y = list_count

            plt.plot(x, y)

            # Save plot
            plt.show()
            plt.savefig(name_fig)
            plt.close()
    except Exception as e:
        print('BAM_file: ', BAM_file)
        print('Exception:', e)
    
    sys.stdout.flush()
    return(BAM_file, datetime.now() - t0)



if __name__=='__main__':
#     if len(sys.argv)<2:
#         print("Usage: {} {}".format(sys.argv[0], "basefolder"))
#         print("Example:")
#         print("       python {} {} {} {}".format(sys.argv[0], "/home/madhu/datasets/download/m6A_RNA_modification_native_RNA_seq/GSM3528751/extracted_data"))
# #       print("       python {} {} {} {}".format(sys.argv[0], "Curlcake_non"))
#         sys.exit(1)

    check_env('ACONDA')

    startTime = datetime.now()
    current_time = startTime.strftime("%Y/%m/%d at %H:%M:%S")
    print('\n\nThe script starting at: ' + str(current_time), ' \n\n' )

    base_folder = sys.argv[1]
    base_folder = check_endWith(base_folder)

    all_BAM_files = [os.path.join(root, name)
                for root, dirs, files in os.walk(base_folder)
                for name in files
                if name.endswith((".bam"))]
    
    inputs = zip(all_BAM_files)
    # pbar = tqdm(total=len(all_BAM_files))

    ## Use the maximum CPU possible 
    cpus = cpu_count()
    n_jobs = min(cpus, len(all_BAM_files))

    # n_jobs = 20

    print('n_jobs: ', n_jobs)
    ## I am using multiple Processors as it is faster than threads. 
    results = tqdm(Pool(n_jobs).imap_unordered(do_MAPQ, inputs), total=len(all_BAM_files))
    # results = ThreadPool(n_jobs).imap_unordered(unzip, inputs)
    for result in results:
        # print('time (s):', result)
        print('File: ', result[0], 'time :', result[1])



    # for BAM_file in tqdm(all_BAM_files[0:1]):
        # output_dir = BAM_file.replace(BAM_file.split('/')[ len(BAM_file.split('/'))-1], '')
        # # output_dir = BAM_file.replace(BAM_file.split('/')[ len(BAM_file.split('/'))-1], 'long_output/') 
        # print('BAM_file: ', BAM_file)
        # print('out_dir: ', output_dir)
        # print('\n')
        # # if check_create_dir(output_dir):
        # #     longRD_cmd = 'python /home/madhu/datasets/codes/LongReadSum bam -i {} -o {}'.format(BAM_file, output_dir)
        # #     # print('longRD_cmd: ', longRD_cmd)
        # #     os.system(longRD_cmd)
        # # else:
        # #     print('Analysis has already been completed for: ', output_dir)
        # # file = sys.argv[1]
        
        # samtool_cmd = 'samtools view {} | cut -f5 | sort -n | uniq -c'.format(BAM_file)
        # stream = os.popen(samtool_cmd)
        # output = stream.read()
        # print('samtool_cmd is: ', samtool_cmd)
        # print('\n\nThe MAPQ distribution for the file: ', BAM_file, ' is: \n', (output))

        # # print(output[0:4])

        # # print(output.split('\n')[0].split(' '))
        # # print(output.split('\n')[1].split(' '))
        # # print(output.split('\n')[2].split(' '))
        # # print(output.split('\n')[len(output.split('\n'))-2].split(' '))

        # ## The outputs are: 

        # '''
        # (ACONDA) madhu@qgenlabgpu:~$ python /home/madhu/datasets/scripts/chris_scripts/check_MAPQ.py
        # ['', '', '14629', '0']
        # ['', '', '', '', '110', '1']
        # ['', '', '', '', '', '46', '10']
        # ['', '', '', '', '', '32', '9']

        # '''
        
        
        
        # list_count = []
        # list_mapQ = []

        # for a_line in output.split('\n'):
        #     if a_line == '':
        #         continue
        #     s_line = a_line.split(' ')
        #     while('' in s_line):
        #         s_line.remove('')
        # #     print('line: ', s_line)
        #     count = s_line[0]
        #     mapQ = s_line[1]
        # #     print('count: ', count)
        # #     print('mapQ: ', mapQ)
        # #     print('\n\n')
        #     list_count.append(int(count))
        #     list_mapQ.append(int(mapQ))

        # ## Now, as I have sorted out the array in the terminal, no need to sort out in the Python ... 
        # x = list_mapQ
        # y = list_count

        # plt.plot(x, y)

        # # Show plot
        # plt.show()

        # # name = file.replace('.bam', '.png')

        # name = BAM_file.split('.bam')[0]+'_MAPQ.png'

        # plt.savefig(name)
        # plt.close()
    
    executionTime = (datetime.now() - startTime)
    current_time = datetime.now().strftime("%Y/%m/%d at %H:%M:%S")
    print('\n\nThe script completed at: ' + str(current_time))
    print('Execution time: ' + str(executionTime), ' \n\n')



