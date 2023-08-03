
import os,sys
# import time, datetime
sys.path.insert(0, '/home/madhu/work/codes')
from main import check_endWith, check_create_dir, check_env, read_file, write_file, count_lines, save_file, load_file

# check_env('torch')

# import TS_finetune_pred_MP
from datetime import datetime
import glob

import warnings
# warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore")

from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm
import pickle

from natsort import natsorted
from multiprocessing.pool import ThreadPool, Pool
from multiprocessing import cpu_count


'''
This script .... 



'''

def get_scores(pred_output_dir):
    y_score_Methyl = []
    y_score_noMethyl = []
    pred_results = load_file(pred_output_dir+'pred_results.pkl')

    for mapc in pred_results:
        # print("Total sites for <{}>: +{} -{}\n".format( mapc, (len(pred_results[mapc]['+']) if '+' in pred_results[mapc] else 0), (len(pred_results[mapc]['-']) if '-' in pred_results[mapc] else 0) ) );
        for mapstn in pred_results[mapc]:
            for ref_pos in pred_results[mapc][mapstn]:
                
                t_pred = pred_results[mapc][mapstn][ ref_pos ]
                # print("Pred: {}:{}:{} {} {} {} {}".format( mapc, mapstn, ref_pos, t_pred[0][0], t_pred[0][1], t_pred[1][0], t_pred[1][1] ) )
                conf_mat = [t_pred[0][0], t_pred[0][1], t_pred[1][0], t_pred[1][1]]
                if (int(conf_mat[0])+int(conf_mat[1])) != 0 and (int(conf_mat[2])+int(conf_mat[3])) != 0:
                    # y_true_ours.append(0)
                    y_score_noMethyl.append(int(conf_mat[1])/(int(conf_mat[0])+int(conf_mat[1])))
                    # y_true_ours.append(1)
                    y_score_Methyl.append(int(conf_mat[3])/(int(conf_mat[2])+int(conf_mat[3])))
    return y_score_Methyl, y_score_noMethyl

def do_get_READ_scores(args):
    # t0 = datetime.now()
    a_part = args[0]
    Nmethyl_read_scores_part = []
    methyl_read_scores_part = []
    if a_part != '':
        lines = a_part.split('\n')
        # print((lines))
        
        pos_list = []
        label_list = []
        actual_label_list = []
        actual_label = 'X'
        for i in range(len(lines)):
            # print(i)
            if lines[i] != '':
                
                if i == 0:
                    # print(lines[i])
                    actual_label = int(lines[i].split(' ')[0])
                    actual_label_list.append(actual_label)
                    readID = lines[i].split(' ')[1]
                    # print('readID: ', readID)
                    # read_list.append(readID)
                else:
                    pos = int(lines[i].split(' ')[0])
                    # print('pos: ', pos)
                    label = int(lines[i].split(' ')[1])
                    # print('label: ', label)
                    pos_list.append(pos)
                    label_list.append(label)
        if len(label_list) != 0:
            # print('label_list: ', label_list)
            # read_scores.append(label_list.count(1)/len(label_list))
            if actual_label == 0:
                # Nmethyl_read_list.append(readID)
                Nmethyl_read_scores_part.append(label_list.count(1)/len(label_list))
            elif actual_label == 1:
                # methyl_read_list.append(readID)
                methyl_read_scores_part.append(label_list.count(1)/len(label_list))


    # print('Nmethyl_read_list: ', Nmethyl_read_list)
    # print('Nmethyl_read_scores: ', Nmethyl_read_scores)

    # print('methyl_read_list: ', methyl_read_list)
    # print('methyl_read_scores: ', methyl_read_scores)

    return methyl_read_scores_part, Nmethyl_read_scores_part



def get_READ_scores(pred_output_dir):
        ## read the predres file ... 
    predres_content = read_file(glob.glob(pred_output_dir+'*.predres')[0])
    # print(predres_content.split('>')[1].split(' ')[1])
    # print(predres_content.split('>')[1].split('\n')[0].split(' ')[1])
    # print(predres_content.split('>')[1].split('\n')[1])
    read_list = []
    Nmethyl_read_list = []
    methyl_read_list = []
    Nmethyl_read_scores = []
    methyl_read_scores = []

    list_sections = predres_content.split('>')

    inputs = zip(list_sections)

    ## Use the maximum CPU possible 
    cpus = cpu_count()
    # n_jobs = min(cpus, len(list_sections))
    n_jobs = 20

    print('n_jobs: ', n_jobs)
    ## I am using multiple Processors as it is faster than threads. 
    results = tqdm(Pool(n_jobs).imap_unordered(do_get_READ_scores, inputs), total=len(list_sections))
    # results = Pool(n_jobs).imap_unordered(do_get_READ_scores, inputs)
    # results = ThreadPool(n_jobs).imap_unordered(unzip, inputs)
    
    
    for result in results:
        methyl_read_scores.extend(result[0])
        Nmethyl_read_scores.extend(result[1])
        # print('time (s):', result)
        # print('methyl_read_scores: ', result[0], 'Nmethyl_read_scores :', result[1])
    


    # # for a_part in tqdm(predres_content.split('>')[0:5]):
    # for a_part in tqdm(predres_content.split('>')):
    #     if a_part != '':
    #         lines = a_part.split('\n')
    #         # print((lines))
            
    #         pos_list = []
    #         label_list = []
    #         actual_label_list = []
    #         actual_label = 'X'
    #         for i in range(len(lines)):
    #             # print(i)
    #             if lines[i] != '':
                    
    #                 if i == 0:
    #                     # print(lines[i])
    #                     actual_label = int(lines[i].split(' ')[0])
    #                     actual_label_list.append(actual_label)
    #                     readID = lines[i].split(' ')[1]
    #                     # print('readID: ', readID)
    #                     # read_list.append(readID)
    #                 else:
    #                     pos = int(lines[i].split(' ')[0])
    #                     # print('pos: ', pos)
    #                     label = int(lines[i].split(' ')[1])
    #                     # print('label: ', label)
    #                     pos_list.append(pos)
    #                     label_list.append(label)
    #         if len(label_list) != 0:
    #             # print('label_list: ', label_list)
    #             # read_scores.append(label_list.count(1)/len(label_list))
    #             if actual_label == 0:
    #                 Nmethyl_read_list.append(readID)
    #                 Nmethyl_read_scores.append(label_list.count(1)/len(label_list))
    #             elif actual_label == 1:
    #                 methyl_read_list.append(readID)
    #                 methyl_read_scores.append(label_list.count(1)/len(label_list))


    # # print('Nmethyl_read_list: ', Nmethyl_read_list)
    # # print('Nmethyl_read_scores: ', Nmethyl_read_scores)

    # # print('methyl_read_list: ', methyl_read_list)
    # # print('methyl_read_scores: ', methyl_read_scores)

    return methyl_read_scores, Nmethyl_read_scores

def do_plot_pred(pred_output_dir_list, out_file, pred_type):
    fig, ax = plt.subplots(1, 3, figsize=FIG_SIZE, sharey=True)
    for i in range(len(pred_output_dir_list)):
        if pred_type == 'SITE':
            y_score_Methyl, y_score_noMethyl = get_scores(pred_output_dir_list[i])
        elif pred_type == 'READ':
            y_score_Methyl, y_score_noMethyl = get_READ_scores(pred_output_dir_list[i])
        hist_KO,bins_KO=np.histogram(np.array(y_score_Methyl),bins=np.linspace(0,1,N_BIN))
        hist_WT,bins_WT=np.histogram(np.array(y_score_noMethyl),bins=np.linspace(0,1,N_BIN))
        print('hist_KO: ', hist_KO)
        print('hist_WT: ', hist_WT)
        # print('sum(hist): ', sum(hist))
        print('bins_KO: ', bins_KO)
        print('bins_WT: ', bins_WT)
        
        # plt.clf()
        ax[i].hist([ np.array(y_score_Methyl) , np.array(y_score_noMethyl) ], bins=np.linspace(0,1,N_BIN), histtype='bar', label=['Methyl Data', 'noMethyl Data'])
        ax[i].legend(loc='best', fontsize=18)
        ax[i].set(xticks=[x for x in np.linspace(0,1,N_BIN)], xticklabels=[round(x, 2) for x in np.linspace(0,1,N_BIN)])
        # ax[0].set_xticklabels([round(x, 2) for x in np.linspace(0,1,N_BIN)], fontsize=20)
        # ax[0].set_yticks(fontsize=20)
        ax[i].grid(color='y', linestyle='-', linewidth=1)
        ax[i].set_xlabel('Methylation Probablity Range', fontsize=22)
        ax[i].set_ylabel('Number of sites', fontsize=22)
        ax[i].set_title('iteration: '+str(pred_output_dir_list[i].split('.pt__')[0].split('.')[-1]), fontsize=25)

    plt.show()
    plt.savefig(out_file)
    plt.close()
    print('\nThe figure has been successfully generated as: ', out_file, '\n')


def plot_pred_OLD():
    ## Plot SITE based scores
    fig, ax = plt.subplots(1, 3, figsize=FIG_SIZE, sharey=True)
    for i in range(len(pred_output_dir_list)):
        y_score_Methyl, y_score_noMethyl = get_scores(pred_output_dir_list[i])
        
        hist_KO,bins_KO=np.histogram(np.array(y_score_Methyl),bins=np.linspace(0,1,N_BIN))
        hist_WT,bins_WT=np.histogram(np.array(y_score_noMethyl),bins=np.linspace(0,1,N_BIN))
        print('hist_KO: ', hist_KO)
        print('hist_WT: ', hist_WT)
        # print('sum(hist): ', sum(hist))
        print('bins_KO: ', bins_KO)
        print('bins_WT: ', bins_WT)
        
        # plt.clf()
        ax[i].hist([ np.array(y_score_Methyl) , np.array(y_score_noMethyl) ], bins=np.linspace(0,1,N_BIN), histtype='bar', label=['Methyl Data', 'noMethyl Data'])
        ax[i].legend(loc='best', fontsize=18)
        ax[i].set(xticks=[x for x in np.linspace(0,1,N_BIN)], xticklabels=[round(x, 2) for x in np.linspace(0,1,N_BIN)])
        # ax[0].set_xticklabels([round(x, 2) for x in np.linspace(0,1,N_BIN)], fontsize=20)
        # ax[0].set_yticks(fontsize=20)
        ax[i].grid(color='y', linestyle='-', linewidth=1)
        ax[i].set_xlabel('Methylation Probablity Range', fontsize=22)
        ax[i].set_ylabel('Number of sites', fontsize=22)
        ax[i].set_title('iteration: '+str(pred_output_dir_list[i].split('.pt__')[0].split('.')[-1]), fontsize=25)

    plt.show()
    plt.savefig(output_PNG_filename)
    plt.close()
    print('\nThe figure has been successfully generated as: ', output_PNG_filename, '\n')
    







    ## Plot READ based scores
    fig, ax = plt.subplots(1, len(pred_output_dir_list), figsize=FIG_SIZE, sharey=True)
    for i in range(len(pred_output_dir_list)):
        y_score_Methyl, y_score_noMethyl = get_READ_scores(pred_output_dir_list[i])
        hist_KO,bins_KO=np.histogram(np.array(y_score_Methyl),bins=np.linspace(0,1,N_BIN))
        hist_WT,bins_WT=np.histogram(np.array(y_score_noMethyl),bins=np.linspace(0,1,N_BIN))
        print('hist_KO: ', hist_KO)
        print('hist_WT: ', hist_WT)
        # print('sum(hist): ', sum(hist))
        print('bins_KO: ', bins_KO)
        print('bins_WT: ', bins_WT)
        
        # plt.clf()
        ax[i].hist([ np.array(y_score_Methyl) , np.array(y_score_noMethyl) ], bins=np.linspace(0,1,N_BIN), histtype='bar', label=['Methyl Data', 'noMethyl Data'])
        ax[i].legend(loc='best', fontsize=18)
        ax[i].set(xticks=[x for x in np.linspace(0,1,N_BIN)], xticklabels=[round(x, 2) for x in np.linspace(0,1,N_BIN)])
        # ax[0].set_xticklabels([round(x, 2) for x in np.linspace(0,1,N_BIN)], fontsize=20)
        # ax[0].set_yticks(fontsize=20)
        ax[i].grid(color='y', linestyle='-', linewidth=1)
        ax[i].set_xlabel('Methylation Probablity Range', fontsize=22)
        ax[i].set_ylabel('Number of Reads', fontsize=22)

    plt.show()
    plt.savefig(output_READ_PNG_filename)
    plt.close()
    print('\nThe figure has been successfully generated as: ', output_READ_PNG_filename, '\n')



## DO NOT DELETE: this one plots one figure, without subplotting 
# def plot_hist(pred_output_dir):
#     ## This is using the saved pickle file .................................. (for real data: mouse)
    
#     y_score_ours_KO, y_score_ours_WT = get_scores(pred_output_dir)

#     # print('The length of sites: ', len(y_score_ours))
#     N_BIN = 20
#     hist_KO,bins_KO=np.histogram(np.array(y_score_ours_KO),bins=np.linspace(0,1,N_BIN))
#     hist_WT,bins_WT=np.histogram(np.array(y_score_ours_WT),bins=np.linspace(0,1,N_BIN))
#     print('hist_KO: ', hist_KO)
#     print('hist_WT: ', hist_WT)
#     # print('sum(hist): ', sum(hist))
#     print('bins_KO: ', bins_KO)
#     print('bins_WT: ', bins_WT)

#     fig, ax = plt.subplots(figsize=(20, 20))
#     plt.clf()
#     # plt.hist(np.array(y_score_ours), bins=np.linspace(0,1,N_BIN), histtype='bar', label='y_score distribution')
#     plt.hist([ np.array(y_score_ours_KO) , np.array(y_score_ours_WT) ], bins=np.linspace(0,1,N_BIN), histtype='bar', label=['KO', 'WT'])
#     plt.legend(loc='best', fontsize=25)
#     plt.xticks([x for x in np.linspace(0,1,N_BIN)], [round(x, 2) for x in np.linspace(0,1,N_BIN)], fontsize=20)
#     plt.yticks(fontsize=20)
#     plt.grid(color='y', linestyle='-', linewidth=1)
#     plt.xlabel('Probablity Thresholds', fontsize=22)
#     plt.ylabel('Number of sites', fontsize=22)
#     plt.show()
    
#     output_PNG_filename = pred_output_dir+ 'pred_SITE_dist.png'
#     plt.savefig(output_PNG_filename)
#     plt.close()
#     print('\nThe figure has been successfully generated as: ', output_PNG_filename, '\n')

def pred_plot():
    pred_output_dir_list = [pred_output_dir_1, pred_output_dir_2, pred_output_dir_3]
    do_plot_pred(pred_output_dir_list, output_PNG_filename, 'SITE')
    do_plot_pred(pred_output_dir_list, output_READ_PNG_filename, 'READ')


if __name__=='__main__' and '__file__' in globals():

    startTime = datetime.now()
    current_time = startTime.strftime("%Y/%m/%d at %H:%M:%S")
    print('\n\nThe script starting at: ' + str(current_time), ' \n\n' )

    FIG_SIZE = (50, 10)
    N_BIN = 20


    '''
    Training data: paper 5 IVT m6A 
    '''
    ## Testing on paper 3 
    # output_PNG_filename = '/home/madhu/work/class_outputs/p5_decoding_epitranscriptional_native_RNA_seq/'+ 'compare_pred_SITE_dist.png'
    output_PNG_filename = '/home/madhu/work/class_outputs/p5_decoding_epitranscriptional_native_RNA_seq/'+ 'pred_SITE_dist_train_p5_test_p3.png'
    output_READ_PNG_filename = '/home/madhu/work/class_outputs/p5_decoding_epitranscriptional_native_RNA_seq/'+ 'pred_READ_dist_train_p5_test_p3.png'

    pred_output_dir_1 = '/home/madhu/work/class_outputs/p5_decoding_epitranscriptional_native_RNA_seq/model_PRED_output__pred_p3_train_p5_m6A_ALL_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.90441.pt__A/'
    pred_output_dir_2 = '/home/madhu/work/class_outputs/p5_decoding_epitranscriptional_native_RNA_seq/model_PRED_output__pred_p3_train_p5_m6A_ALL_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.210432.pt__A/'
    pred_output_dir_3 = '/home/madhu/work/class_outputs/p5_decoding_epitranscriptional_native_RNA_seq/model_PRED_output__pred_p3_train_p5_m6A_ALL_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.330287.pt__A/'
    
    pred_plot()


    ## Testing on the same IVT data, used for training 
    # output_PNG_filename = '/home/madhu/work/class_outputs/p5_decoding_epitranscriptional_native_RNA_seq/'+ 'compare_SAME_IVT_pred_SITE_dist.png'
    output_PNG_filename = '/home/madhu/work/class_outputs/p5_decoding_epitranscriptional_native_RNA_seq/'+ 'pred_SITE_dist_train_p5_test_p5.png'
    output_READ_PNG_filename = '/home/madhu/work/class_outputs/p5_decoding_epitranscriptional_native_RNA_seq/'+ 'pred_READ_dist_train_p5_test_p5.png'

    pred_output_dir_1 = '/home/madhu/work/class_outputs/p5_decoding_epitranscriptional_native_RNA_seq/model_PRED_output__pred_p5_train_p5_m6A_overfit_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.90441.pt__A/'
    pred_output_dir_2 = '/home/madhu/work/class_outputs/p5_decoding_epitranscriptional_native_RNA_seq/model_PRED_output__pred_p5_train_p5_m6A_overfit_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.210432.pt__A/'
    pred_output_dir_3 = '/home/madhu/work/class_outputs/p5_decoding_epitranscriptional_native_RNA_seq/model_PRED_output__pred_p5_train_p5_m6A_overfit_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.330287.pt__A/'

    pred_plot()



    '''
    Training data: paper 3 curlcake
    '''
    ## Testing data: paper 5
    output_PNG_filename = '/home/madhu/work/class_outputs/p5_decoding_epitranscriptional_native_RNA_seq/'+ 'pred_SITE_dist_train_p3_test_p5.png'
    output_READ_PNG_filename = '/home/madhu/work/class_outputs/p5_decoding_epitranscriptional_native_RNA_seq/'+ 'pred_READ_dist_train_p3_test_p5.png'

    pred_output_dir_1 = '/home/madhu/work/class_outputs/p3_m6A_RNA_modification_native_RNA_seq/model_PRED_output__train_p3_BOTH_pred_p5_m6A_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.70055.pt__A/'
    pred_output_dir_2 = '/home/madhu/work/class_outputs/p3_m6A_RNA_modification_native_RNA_seq/model_PRED_output__train_p3_BOTH_pred_p5_m6A_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.240025.pt__A/'
    pred_output_dir_3 = '/home/madhu/work/class_outputs/p3_m6A_RNA_modification_native_RNA_seq/model_PRED_output__train_p3_BOTH_pred_p5_m6A_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.650154.pt__A/'

    pred_plot()


    ## Testing data: paper 3
    output_PNG_filename = '/home/madhu/work/class_outputs/p5_decoding_epitranscriptional_native_RNA_seq/'+ 'pred_SITE_dist_train_p3_test_p3.png'
    output_READ_PNG_filename = '/home/madhu/work/class_outputs/p5_decoding_epitranscriptional_native_RNA_seq/'+ 'pred_READ_dist_train_p3_test_p3.png'

    pred_output_dir_1 = '/home/madhu/work/class_outputs/p3_m6A_RNA_modification_native_RNA_seq/model_PRED_output__train_p3_overfitting_pred_p3_m6A_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.70055.pt__A/'
    pred_output_dir_2 = '/home/madhu/work/class_outputs/p3_m6A_RNA_modification_native_RNA_seq/model_PRED_output__train_p3_overfitting_pred_p3_m6A_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.240025.pt__A/'
    pred_output_dir_3 = '/home/madhu/work/class_outputs/p3_m6A_RNA_modification_native_RNA_seq/model_PRED_output__train_p3_overfitting_pred_p3_m6A_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.650154.pt__A/'

    pred_plot()


    # pred_output_dir_list = [pred_output_dir_1, pred_output_dir_2, pred_output_dir_3]

    # plot_pred(output_PNG_filename, 'SITE')
    # plot_pred(output_READ_PNG_filename, 'READ')


    # ## Plot SITE based scores
    # fig, ax = plt.subplots(1, 3, figsize=FIG_SIZE, sharey=True)
    # for i in range(len(pred_output_dir_list)):
    #     y_score_Methyl, y_score_noMethyl = get_scores(pred_output_dir_list[i])
    #     N_BIN = 20
    #     hist_KO,bins_KO=np.histogram(np.array(y_score_Methyl),bins=np.linspace(0,1,N_BIN))
    #     hist_WT,bins_WT=np.histogram(np.array(y_score_noMethyl),bins=np.linspace(0,1,N_BIN))
    #     print('hist_KO: ', hist_KO)
    #     print('hist_WT: ', hist_WT)
    #     # print('sum(hist): ', sum(hist))
    #     print('bins_KO: ', bins_KO)
    #     print('bins_WT: ', bins_WT)
        
    #     # plt.clf()
    #     ax[i].hist([ np.array(y_score_Methyl) , np.array(y_score_noMethyl) ], bins=np.linspace(0,1,N_BIN), histtype='bar', label=['Methyl Data', 'noMethyl Data'])
    #     ax[i].legend(loc='best', fontsize=18)
    #     ax[i].set(xticks=[x for x in np.linspace(0,1,N_BIN)], xticklabels=[round(x, 2) for x in np.linspace(0,1,N_BIN)])
    #     # ax[0].set_xticklabels([round(x, 2) for x in np.linspace(0,1,N_BIN)], fontsize=20)
    #     # ax[0].set_yticks(fontsize=20)
    #     ax[i].grid(color='y', linestyle='-', linewidth=1)
    #     ax[i].set_xlabel('Methylation Probablity Range', fontsize=22)
    #     ax[i].set_ylabel('Number of sites', fontsize=22)
    #     ax[i].set_title('iteration: '+str(pred_output_dir_list[i].split('.pt__')[0].split('.')[-1]), fontsize=25)

    # plt.show()
    # plt.savefig(output_PNG_filename)
    # plt.close()
    # print('\nThe figure has been successfully generated as: ', output_PNG_filename, '\n')
    







    # ## Plot READ based scores
    # fig, ax = plt.subplots(1, len(pred_output_dir_list), figsize=FIG_SIZE, sharey=True)
    # for i in range(len(pred_output_dir_list)):
    #     y_score_Methyl, y_score_noMethyl = get_READ_scores(pred_output_dir_list[i])
    #     N_BIN = 20
    #     hist_KO,bins_KO=np.histogram(np.array(y_score_Methyl),bins=np.linspace(0,1,N_BIN))
    #     hist_WT,bins_WT=np.histogram(np.array(y_score_noMethyl),bins=np.linspace(0,1,N_BIN))
    #     print('hist_KO: ', hist_KO)
    #     print('hist_WT: ', hist_WT)
    #     # print('sum(hist): ', sum(hist))
    #     print('bins_KO: ', bins_KO)
    #     print('bins_WT: ', bins_WT)
        
    #     # plt.clf()
    #     ax[i].hist([ np.array(y_score_Methyl) , np.array(y_score_noMethyl) ], bins=np.linspace(0,1,N_BIN), histtype='bar', label=['Methyl Data', 'noMethyl Data'])
    #     ax[i].legend(loc='best', fontsize=18)
    #     ax[i].set(xticks=[x for x in np.linspace(0,1,N_BIN)], xticklabels=[round(x, 2) for x in np.linspace(0,1,N_BIN)])
    #     # ax[0].set_xticklabels([round(x, 2) for x in np.linspace(0,1,N_BIN)], fontsize=20)
    #     # ax[0].set_yticks(fontsize=20)
    #     ax[i].grid(color='y', linestyle='-', linewidth=1)
    #     ax[i].set_xlabel('Methylation Probablity Range', fontsize=22)
    #     ax[i].set_ylabel('Number of Reads', fontsize=22)

    # plt.show()
    # plt.savefig(output_READ_PNG_filename)
    # plt.close()
    # print('\nThe figure has been successfully generated as: ', output_READ_PNG_filename, '\n')



    
    executionTime = (datetime.now() - startTime)
    current_time = datetime.now().strftime("%Y/%m/%d at %H:%M:%S")
    print('\n\nThe script completed at: ' + str(current_time))
    print('Execution time: ' + str(executionTime), ' \n\n')

