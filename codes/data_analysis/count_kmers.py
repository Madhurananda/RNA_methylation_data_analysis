
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 11:35:53 2023

@author: pahar
"""

'''
This script analyses the spike protiens ... 
'''

import os
import sys
sys.path.insert(0, '/home/madhu/work/codes')

from main import check_endWith, check_create_dir, check_env, del_dir_content_files, read_file

import glob
# import librosa
import pandas as pd
# import soundfile as sf
import matplotlib.pyplot as plt
import numpy as np
# from natsort import natsorted
from tqdm import tqdm 
from datetime import datetime
import itertools




############################
# Download the data and load it into Python 

'''

Data has been downloaded from the following locations: 

Dataset 1: 
https://www.ncbi.nlm.nih.gov/gene/1489668

Dataset 2:
https://www.ncbi.nlm.nih.gov/gene/43740568

Dataset 3:
https://www.ncbi.nlm.nih.gov/gene/14254594/

'''

#########################
# There are a few types of data and I am looking into the .fna and .faa 
#########################


# data_dir = '../data/S_datasets_3/ncbi_dataset/data/'



# faa_file = list(glob.iglob(data_dir+'*.faa*', recursive=True))

# fna_file = list(glob.iglob(data_dir+'*.fna*', recursive=True))



# f = open(faa_file[0], "r")
# # print(f.read())
# faa_content = f.read()
# f.close()

# f = open(fna_file[0], "r")
# fna_content = f.read()
# # print(f.read())
# f.close()



### Now it seems that I need to use the FNA file ... 



def count_kmers(read, k):
    """Count kmer occurrences in a given read.

    Parameters
    ----------
    read : string
        A single DNA sequence.
    k : int
        The value of k for which to count kmers.

    Returns
    -------
    counts : dictionary, {'string': int}
        A dictionary of counts keyed by their individual kmers (strings
        of length k).

    Examples
    --------
    >>> count_kmers("GATGAT", 3)
    {'ATG': 1, 'GAT': 2, 'TGA': 1}
    """
    # Start with an empty dictionary
    counts = {}
    # Calculate how many kmers of length k there are
    num_kmers = len(read) - k + 1
    # Loop over the kmer start positions
    for i in range(num_kmers):
        # Slice the string to get the kmer
        kmer = read[i:i+k]
        # Add the kmer to the dictionary if it's not there
        if kmer not in counts:
            counts[kmer] = 0
        # Increment the count for this kmer
        counts[kmer] += 1
    # Return the final counts
    return counts




if __name__=='__main__':
#     if len(sys.argv)<2:
#         print("Usage: {} {}".format(sys.argv[0], "basefolder"))
#         print("Example:")
#         print("       python {} {} {} {}".format(sys.argv[0], "/home/madhu/datasets/download/m6A_RNA_modification_native_RNA_seq/GSM3528751/extracted_data"))
# #           print("       python {} {} {} {}".format(sys.argv[0], "Curlcake_non"))
#         sys.exit(1)
    
    check_env('ACONDA')

    startTime = datetime.now()
    current_time = startTime.strftime("%Y/%m/%d at %H:%M:%S")
    print('\n\nThe script starting at: ' + str(current_time), ' \n\n' )


    input_file = sys.argv[1]
    file_type = sys.argv[2]

    print('Counting 5-mers ... ')
    if file_type == 'REF':
        print('')
        file_content = ''
        # fig_c = 1

        for a_line in tqdm(read_file(input_file).split('\n')):
            if a_line.startswith('A') or a_line.startswith('C') or a_line.startswith('G') or a_line.startswith('T'):
                # print('')
                file_content = a_line + file_content
        
        mer_5_count = count_kmers(file_content,5)

    else:
        print('The file types should be: REF or ...')
        sys.exit()



    # print('mer_5_count')
    # print(mer_5_count)




    alphabets = ['A', 'C', 'G', 'T']

    final_combinations = []

    combinations = itertools.product(*itertools.repeat(alphabets, 5))
    for i in combinations:
        # print(i)
        letter = ''
        for j in i:
            letter = j+letter
        final_combinations.append(letter)


    final_dict = {}

    for a_comb in final_combinations:
        # print(a_comb)
        if a_comb in mer_5_count:
            # print('')
            final_dict[a_comb] = mer_5_count[a_comb]
        else:
            final_dict[a_comb] = 0


    # df = pd.DataFrame.from_dict(final_dict)

    df = pd.DataFrame.from_dict(final_dict, orient='index', columns=['Count'])


    # df.sort_values(by=['Count'], inplace=True)

    df.sort_values(by=['Count'], inplace=True, ascending=False)


    # print(df['Count'].value_counts()[3])



    # saveFig_title = '../results/freq_plot'
    saveFig_title = input_file.replace('.'+input_file.split('.')[-1], '_5mer.png')

    fig, ax = plt.subplots(figsize=(20, 10))
    plt.plot(df['Count'].values, 'b.')
    # plt.plot(df['Count'].values)


    plt.xlabel('Nucleotides', fontsize=20)
    plt.ylabel('5-mer Frequency', fontsize=20)
    plt.grid(color='y', linewidth=0.5)
    # plt.xlim([-0.02, 1050])
    # plt.ylim([-0.02, 23])
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.title('Distribution of 5-mers', fontsize=35)
    plt.show()


    plt.savefig(saveFig_title)
    plt.close('all')





    # csv_save_file = '../results/5-mer-results.csv'
    csv_save_file = input_file.replace('.'+input_file.split('.')[-1], '_5mer.csv')

    df.to_csv(csv_save_file, index=True)





    df_5mer_0 = df[df['Count']== 0]

    # csv_save_file = '../results/5-mer-results__0.csv'
    csv_save_file = input_file.replace('.'+input_file.split('.')[-1], '_5mer__0.csv')
    df_5mer_0.to_csv(csv_save_file, index=True)



    df_5mer_1 = df[df['Count']== 1]

    # csv_save_file = '../results/5-mer-results__1.csv'
    csv_save_file = input_file.replace('.'+input_file.split('.')[-1], '_5mer__1.csv')
    df_5mer_1.to_csv(csv_save_file, index=True)








    '''
    Plot the best 10/25 nucleotides
    '''



    top_10_df = df[0:25]

    X_Values = top_10_df.index

    Y_Values = top_10_df['Count'].values




    saveFig_title = csv_save_file = input_file.replace('.'+input_file.split('.')[-1], '_5mer_freq_plot_best_25.png')

    '../results/freq_plot_best_25'

    fig, ax = plt.subplots(figsize=(20, 10))
    plt.stem(X_Values, Y_Values, 'bo')
    # plt.plot(df['Count'].values)


    plt.xlabel('Nucleotides', fontsize=20)
    plt.ylabel('5-mer Frequency', fontsize=20)
    plt.grid(color='y', linewidth=0.5)
    # plt.xlim([-0.02, 1050])
    # plt.ylim([10.5, 20.5])
    plt.ylim([ int(min(Y_Values)*0.95), int(max(Y_Values)*1.05) ])



    # n_split = 11
    # xTick_list = []
    # for n in np.linspace(0, 1, n_split):
    #     xTick_list.append(str(int(n*100))+'%')
    # plt.xticks(np.linspace(0, 1, n_split), xTick_list, fontsize=15)
    # yTick_list = []
    # for n in np.linspace(0, 1, n_split):
    #     yTick_list.append(str(int(n*100))+'%')
    # plt.yticks(np.linspace(0, 1, n_split), yTick_list, fontsize=15)


    plt.xticks(fontsize=15, rotation=30)
    plt.yticks(fontsize=15)
    plt.title('25 Most frequent 5-mers', fontsize=35)
    plt.show()


    plt.savefig(saveFig_title)
    plt.close('all')



    executionTime = (datetime.now() - startTime)
    current_time = datetime.now().strftime("%Y/%m/%d at %H:%M:%S")
    print('\n\nThe script completed at: ' + str(current_time))
    print('Execution time: ' + str(executionTime), ' \n\n')


