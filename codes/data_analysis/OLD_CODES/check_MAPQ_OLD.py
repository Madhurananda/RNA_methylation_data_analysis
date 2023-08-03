
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from matplotlib.ticker import PercentFormatter

from operator import itemgetter


'''
I am creating this file to cehck the MAPQ (mapping quality) in a BAM file. 
This file should produce a distribution of all the MAPQ values. 
'''


if __name__=='__main__':
#     if len(sys.argv)<2:
#         print("Usage: {} {}".format(sys.argv[0], "basefolder"))
#         print("Example:")
#         print("       python {} {} {} {}".format(sys.argv[0], "/home/madhu/datasets/download/m6A_RNA_modification_native_RNA_seq/GSM3528751/extracted_data"))
# #       print("       python {} {} {} {}".format(sys.argv[0], "Curlcake_non"))
#         sys.exit(1)
   
    file = sys.argv[1]

    
#     file="/home/madhu/datasets/analysis/p3_m6A_RNA_modification_native_RNA_seq/GSM3528749/GSM3528749.bam"
    # samtools view "$file" | cut -f5 | sort | uniq -c


    stream = os.popen('samtools view {} | cut -f5 | sort -n | uniq -c'.format(file))
    output = stream.read()

    print('\n\nThe MAPQ distribution for the file: ', file, ' is: \n', (output))

    # print(output[0:4])

    # print(output.split('\n')[0].split(' '))
    # print(output.split('\n')[1].split(' '))
    # print(output.split('\n')[2].split(' '))
    # print(output.split('\n')[len(output.split('\n'))-2].split(' '))

    ## The outputs are: 

    '''
    (ACONDA) madhu@qgenlabgpu:~$ python /home/madhu/datasets/scripts/chris_scripts/check_MAPQ.py
    ['', '', '14629', '0']
    ['', '', '', '', '110', '1']
    ['', '', '', '', '', '46', '10']
    ['', '', '', '', '', '32', '9']

    '''
    
    
    
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


#     ## Just sorting the array in Python .... 
#     a = np.array(list_mapQ).reshape(len(list_count), 1)
#     b = np.array(list_count).reshape(len(list_mapQ), 1)
#     c = np.hstack((a,b))
#     d = sorted(list(map(list, c)), key = itemgetter(0))
    
#     # Creating distribution
#     x = [x[0] for x in d]
#     y = [x[1] for x in d]

    ## Now, as I have sorted out the array in the terminal, no need to sort out in the Python ... 
    x = list_mapQ
    y = list_count

    plt.plot(x, y)

    # Show plot
    plt.show()

    # name = file.replace('.bam', '.png')

    name = file.split('.bam')[0]+'_MAPQ.png'

    plt.savefig(name)
    plt.close()




