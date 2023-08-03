
import os
import sys

import matplotlib.pyplot as plt
import math
import numpy as np
from tqdm import tqdm

'''
This script plots the predres files.  
The basefolder neeeds to be provided and it will scan thourgh all the .predres files and create figures for each of them. 

### Example call: 
python /home/madhu/work/codes/ML_codes/madhu_codes/predres_plots.py /home/madhu/work/class_outputs/p3_m6A_RNA_modification_native_RNA_seq

'''


def draw_plot(predres_file, output_PNG_filename):
    res = open(predres_file,"r")
    Lines = res.readlines()
    res_neg = {}
    res_pos = {}
    res = {}
    p = 0
    pos = 0
    _all = 0
    r = ""
    for line in Lines:
        if len(line) == 1:
            continue
        elif line[0] == ">":
            if line[1] == '0':
                if _all == 0:
                    continue
                p = 100*(pos/_all)
                pos = 0
                _all = 0
                res_neg[r] = p
                r = line.split(" ")[1]
            elif line[1] == '1':
                if _all == 0:
                    continue
                p = 100*(pos/_all)
                pos = 0
                _all = 0
                res_pos[r] = p
                r = line.split(" ")[1]
        else:
            pos = pos + 1.0 if line.strip()[-1] == "0" else pos
            _all += 1.0
    l_neg = []
    for e in res_neg:
        l_neg.append(res_neg[e])
    l_pos = []
    for e in res_pos:
        l_pos.append(res_pos[e])
    

    fig, ax = plt.subplots(figsize=(20, 10))
    plt.clf()
    plt.hist([l_neg, l_pos], bins=[x for x in range(0,101,5)], histtype='bar', label=['Negative', 'Positive'])
    plt.legend(loc='best', fontsize=20)
    plt.xticks([x for x in range(0,101,5)], [x for x in range(0,101,5)], fontsize=15)
    plt.yticks(fontsize=15)
    plt.grid(color='y', linestyle='-', linewidth=1)
    plt.xlabel('Percentage of positives and negatives', fontsize=20)
    plt.ylabel('Number of Reads', fontsize=20)
    
    plt.savefig(output_PNG_filename)
    plt.close()
    print('\nThe figure has been successfully generated as: ', output_PNG_filename, '\n')



if __name__=='__main__':
    if len(sys.argv)<2:
        print("Usage: {} {}".format(sys.argv[0], "basefolder"))
        print("Example:")
        print("       python {} {}".format(sys.argv[0], "/home/madhu/work/class_outputs/p3_m6A_RNA_modification_native_RNA_seq"))
        sys.exit(1)
    
    base_dir = sys.argv[1]

    predres_files = [os.path.join(root, name)
                for root, dirs, files in os.walk(base_dir)
                for name in files
                if name.endswith((".predres"))]
    
    print('\n\n**** There are ', len(predres_files), ' .predres files inside: ', base_dir, '****\n\n')

    for a_predres_file in tqdm(predres_files):
        output_PNG_filename = a_predres_file.replace('.predres', '.png')
        print('\n\nGenerating figures for: ')
        print('a_predres_file: ', a_predres_file)
        print('output_PNG_filename: ', output_PNG_filename, '\n\n')
        draw_plot(a_predres_file, output_PNG_filename)

