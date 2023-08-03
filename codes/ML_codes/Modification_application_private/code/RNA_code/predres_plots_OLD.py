
import os
import sys
sys.path.insert(0, '/home/madhu/work/codes')
from main import check_endWith, check_create_dir, check_env

import matplotlib.pyplot as plt
import math
import numpy as np
from tqdm import tqdm

'''
This script plots the predres. 
The folder: /home/madhu/work/class_outputs/p3_m6A_RNA_modification_native_RNA_seq
The training output folder: model_output__train_p3_3_512_1000_256_1__A
The model test folder:
/home/madhu/work/class_outputs/p3_m6A_RNA_modification_native_RNA_seq/model_PRED_output__train_p3_ALL_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.310124.pt__A/

Files to be read look like: 
less /home/madhu/work/class_outputs/p3_m6A_RNA_modification_native_RNA_seq/model_PRED_output__train_p3_ALL_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.310124.pt__A/ps.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.310124.predres

'''


def draw_plot(predres_file, output_PNG_filename):
    # f = ["TS_finetune_pred.May.04.2023.14.BilstmMean.layers3.hs1024.F.lr500.b256.p1000.GPU.ep1.40029.pt.log", "TS_finetune_pred.May.04.2023.14.BilstmMean.layers3.hs1024.F.lr500.b256.p1000.GPU.ep1.50008.pt.log", "TS_finetune_pred.May.04.2023.14.BilstmMean.layers3.hs1024.F.lr500.b256.p1000.GPU.ep1.60012.pt.log", "TS_finetune_pred.May.04.2023.14.BilstmMean.layers3.hs1024.F.lr500.b256.p1000.GPU.ep1.100011.pt.log", "TS_finetune_pred.May.04.2023.14.BilstmMean.layers3.hs1024.F.lr500.b256.p1000.GPU.ep1.70020.pt.log", "TS_finetune_pred.May.04.2023.14.BilstmMean.layers3.hs1024.F.lr500.b256.p1000.GPU.ep1.10071.pt.log", "TS_finetune_pred.May.04.2023.14.BilstmMean.layers3.hs1024.F.lr500.b256.p1000.GPU.ep1.80069.pt.log", "TS_finetune_pred.May.04.2023.14.BilstmMean.layers3.hs1024.F.lr500.b256.p1000.GPU.ep1.20018.pt.log", "TS_finetune_pred.May.04.2023.14.BilstmMean.layers3.hs1024.F.lr500.b256.p1000.GPU.ep1.90012.pt.log"]

    for x in predres_file:
        res = open(x,"r")
        Lines = res.readlines()
        res_neg = {}
        res_pos = {}
        res = {}
        p = 0
        pos = 0
        _all = 0
        r = ""
        # print(Lines[:5])
        for line in Lines:
            if len(line) == 1:
                continue
            elif line[0] == ">":
                if line[1] == '0':
                # if line[1] == label:
                    if _all == 0:
                        continue
                    p = 100*(pos/_all)
                    pos = 0
                    _all = 0
                    res_neg[r] = p
                    # res[r] = p
                    r = line.split(" ")[1]
                elif line[1] == '1':
                    if _all == 0:
                        continue
                    p = 100*(pos/_all)
                    pos = 0
                    _all = 0
                    res_pos[r] = p
                    # res[r] = p
                    r = line.split(" ")[1]
            else:
                pos = pos + 1.0 if line.strip()[-1] == "0" else pos
                _all += 1.0
        # print(x)
        l_neg = []
        for e in res_neg:
            l_neg.append(res_neg[e])
        l_pos = []
        for e in res_pos:
            l_pos.append(res_pos[e])
        # l = []
        # for e in res:
        #     l.append(res[e])

        # neg_counts, edges, bars = plt.hist(l_neg, bins=[x for x in range(0,101,5)], edgecolor="black", label='Negative')
        # pos_counts, edges, bars = plt.hist(l_pos, bins=[x for x in range(0,101,5)], edgecolor="black", label='Positive')
        plt.close()

        fig, ax = plt.subplots(figsize=(20, 10))

        plt.clf()
        # plt.hist(l_neg, bins=[x for x in range(0,101,5)], edgecolor="black", label='Negative')
        # plt.hist(l_pos, bins=[x for x in range(0,101,5)], edgecolor="black", label='Positive')
        # new_counts, edges, bars = plt.hist([l_neg, l_pos], bins=[x for x in range(0,101,5)], stacked=True, density=True, label=['Negative', 'Positive'])

        plt.hist([l_neg, l_pos], bins=[x for x in range(0,101,5)], histtype='bar', label=['Negative', 'Positive'])

        # hist, edges = np.histogram(
        #     l_neg,
        #     bins=[x for x in range(0,101,5)],
        #     range=(0, bin_width*num_bins),
        #     density=False)

        # print('l_neg: ', l_neg)

        # print('(counts[0]): ', (counts[0]))
        # print('(counts[1]): ', (counts[1]))

        # print('(counts[0][-1]): ', (counts[0][-1]))
        # print('(counts[1][-1]): ', (counts[1][-1]))

        # plt.bar_label(bars, fontsize=12)
        # plt.hist([l_neg, l_pos], bins=[x for x in range(0,101,5)], stacked=True, density=True, edgecolor="black", label=['Negative', 'Positive'])
        plt.legend(loc='best', fontsize=20)

        plt.xticks([x for x in range(0,101,5)], [x for x in range(0,101,5)], fontsize=15)
        

        # max_count = max( max(neg_counts), max(pos_counts))

        # # if max(new_counts[0]) > max(new_counts[1]):
        # if max(neg_counts) > max(pos_counts):
        #     print('This means neg has the higher values')
        #     max_new_counts = max(new_counts[0])*1.05
        #     max_counts = int(max(neg_counts))*1.05
        # else:
        #     print('This means pos has the higher values')
        #     max_new_counts = max(new_counts[1])*1.05
        #     max_counts = int(max(pos_counts))*1.05

        # print('counts: ', len(counts))
        # print('new_counts[0]: ', len(new_counts[0]))
        # plt.yticks(np.arange(0,(max(counts)), (max(counts)/10)), np.arange(0,(max(counts)),(max(counts)/10)), fontsize=15)
        # print('max_new_counts: ', max_new_counts)

        # plt.ylim([0, (max_new_counts*1.05)])

        # ytick_locs = np.arange(0, max_new_counts, (max_new_counts/10)) 
        # ytick_labels = np.arange(0, max_counts, (max_counts/10))

        # print('ytick_locs: ', ytick_locs)
        # print('ytick_labels: ', ytick_labels)

        # plt.yticks(ytick_locs, [int(x) for x in ytick_labels], fontsize=15)

        plt.yticks(fontsize=15)
        plt.grid(color='y', linestyle='-', linewidth=1)
        plt.xlabel('Percentage of positives and negatives', fontsize=20)
        plt.ylabel('Number of Reads', fontsize=20)

        # plt.hist(l, bins=[x for x in range(0,101,5)], edgecolor="black")
        plt.savefig(output_PNG_filename)
        print('\nThe figure has been successfully generated: ', output_PNG_filename)



if __name__=='__main__':
    if len(sys.argv)<2:
        print("Usage: {} {}".format(sys.argv[0], "basefolder"))
        print("Example:")
        print("       python {} {} {} {}".format(sys.argv[0], "/home/madhu/work/download/m6A_RNA_modification_native_RNA_seq/GSM3528751/extracted_data"))
#           print("       python {} {} {} {}".format(sys.argv[0], "Curlcake_non"))
        sys.exit(1)
    
    # check_env('torch')
    base_dir = sys.argv[1]

    predres_files = [os.path.join(root, name)
                for root, dirs, files in os.walk(base_dir)
                for name in files
                if name.endswith((".predres"))]
    
    # if len(predres_file)>1:
    #     print('There is something wrong. Multiple predres files inside: ', base_dir)
    #     sys.exit(0)
    for a_predres_file in tqdm(predres_files):
        output_PNG_filename = a_predres_file.replace('.predres', '.png')
        # label = '1'
        draw_plot(a_predres_file, output_PNG_filename)

