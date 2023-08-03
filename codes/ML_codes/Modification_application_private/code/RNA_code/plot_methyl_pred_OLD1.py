
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
                    y_score_Methyl.append(int(conf_mat[1])/(int(conf_mat[0])+int(conf_mat[1])))
                    # y_true_ours.append(1)
                    y_score_noMethyl.append(int(conf_mat[3])/(int(conf_mat[2])+int(conf_mat[3])))
    return y_score_Methyl, y_score_noMethyl


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



if __name__=='__main__' and '__file__' in globals():

    startTime = datetime.now()
    current_time = startTime.strftime("%Y/%m/%d at %H:%M:%S")
    print('\n\nThe script starting at: ' + str(current_time), ' \n\n' )

    '''
    Training data: paper 5 IVT m6A 
    '''
    # ## Testing on paper 3 
    # # output_PNG_filename = '/home/madhu/work/class_outputs/p5_decoding_epitranscriptional_native_RNA_seq/'+ 'compare_pred_SITE_dist.png'
    # output_PNG_filename = '/home/madhu/work/class_outputs/p5_decoding_epitranscriptional_native_RNA_seq/'+ 'pred_SITE_dist_train_p5_test_p3.png'

    # pred_output_dir_1 = '/home/madhu/work/class_outputs/p5_decoding_epitranscriptional_native_RNA_seq/model_PRED_output__pred_p3_train_p5_m6A_ALL_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.90441.pt__A/'
    # pred_output_dir_2 = '/home/madhu/work/class_outputs/p5_decoding_epitranscriptional_native_RNA_seq/model_PRED_output__pred_p3_train_p5_m6A_ALL_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.210432.pt__A/'
    # pred_output_dir_3 = '/home/madhu/work/class_outputs/p5_decoding_epitranscriptional_native_RNA_seq/model_PRED_output__pred_p3_train_p5_m6A_ALL_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.330287.pt__A/'
    
    

    # ## Testing on the same IVT data, used for training 
    # # output_PNG_filename = '/home/madhu/work/class_outputs/p5_decoding_epitranscriptional_native_RNA_seq/'+ 'compare_SAME_IVT_pred_SITE_dist.png'
    # output_PNG_filename = '/home/madhu/work/class_outputs/p5_decoding_epitranscriptional_native_RNA_seq/'+ 'pred_SITE_dist_train_p5_test_p5.png'

    # pred_output_dir_1 = '/home/madhu/work/class_outputs/p5_decoding_epitranscriptional_native_RNA_seq/model_PRED_output__pred_p5_train_p5_m6A_overfit_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.90441.pt__A/'
    # pred_output_dir_2 = '/home/madhu/work/class_outputs/p5_decoding_epitranscriptional_native_RNA_seq/model_PRED_output__pred_p5_train_p5_m6A_overfit_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.210432.pt__A/'
    # pred_output_dir_3 = '/home/madhu/work/class_outputs/p5_decoding_epitranscriptional_native_RNA_seq/model_PRED_output__pred_p5_train_p5_m6A_overfit_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.330287.pt__A/'



    '''
    Training data: paper 3 curlcake
    '''
    # ## Testing data: paper 5
    # output_PNG_filename = '/home/madhu/work/class_outputs/p5_decoding_epitranscriptional_native_RNA_seq/'+ 'pred_SITE_dist_train_p3_test_p5.png'

    # pred_output_dir_1 = '/home/madhu/work/class_outputs/p3_m6A_RNA_modification_native_RNA_seq/model_PRED_output__train_p3_BOTH_pred_p5_m6A_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.70055.pt__A/'
    # pred_output_dir_2 = '/home/madhu/work/class_outputs/p3_m6A_RNA_modification_native_RNA_seq/model_PRED_output__train_p3_BOTH_pred_p5_m6A_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.240025.pt__A/'
    # pred_output_dir_3 = '/home/madhu/work/class_outputs/p3_m6A_RNA_modification_native_RNA_seq/model_PRED_output__train_p3_BOTH_pred_p5_m6A_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.650154.pt__A/'



    ## Testing data: paper 3
    output_PNG_filename = '/home/madhu/work/class_outputs/p5_decoding_epitranscriptional_native_RNA_seq/'+ 'pred_SITE_dist_train_p3_test_p3.png'

    pred_output_dir_1 = '/home/madhu/work/class_outputs/p3_m6A_RNA_modification_native_RNA_seq/model_PRED_output__train_p3_overfitting_pred_p3_m6A_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.70055.pt__A/'
    pred_output_dir_2 = '/home/madhu/work/class_outputs/p3_m6A_RNA_modification_native_RNA_seq/model_PRED_output__train_p3_overfitting_pred_p3_m6A_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.240025.pt__A/'
    pred_output_dir_3 = '/home/madhu/work/class_outputs/p3_m6A_RNA_modification_native_RNA_seq/model_PRED_output__train_p3_overfitting_pred_p3_m6A_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.650154.pt__A/'



    pred_output_dir_list = [pred_output_dir_1, pred_output_dir_2, pred_output_dir_3]
    fig, ax = plt.subplots(1, 3, figsize=(20, 10), sharey=True)

    for i in range(len(pred_output_dir_list)):
        y_score_Methyl, y_score_noMethyl = get_scores(pred_output_dir_list[i])
        N_BIN = 20
        hist_KO,bins_KO=np.histogram(np.array(y_score_Methyl),bins=np.linspace(0,1,N_BIN))
        hist_WT,bins_WT=np.histogram(np.array(y_score_noMethyl),bins=np.linspace(0,1,N_BIN))
        print('hist_KO: ', hist_KO)
        print('hist_WT: ', hist_WT)
        # print('sum(hist): ', sum(hist))
        print('bins_KO: ', bins_KO)
        print('bins_WT: ', bins_WT)
        
        # plt.clf()
        ax[i].hist([ np.array(y_score_Methyl) , np.array(y_score_noMethyl) ], bins=np.linspace(0,1,N_BIN), histtype='bar', label=['Methyl', 'noMethyl'])
        ax[i].legend(loc='best', fontsize=25)
        ax[i].set(xticks=[x for x in np.linspace(0,1,N_BIN)], xticklabels=[round(x, 2) for x in np.linspace(0,1,N_BIN)])
        # ax[0].set_xticklabels([round(x, 2) for x in np.linspace(0,1,N_BIN)], fontsize=20)
        # ax[0].set_yticks(fontsize=20)
        ax[i].grid(color='y', linestyle='-', linewidth=1)
        ax[i].set_xlabel('Probablity Thresholds', fontsize=22)
        ax[i].set_ylabel('Number of sites', fontsize=22)

    plt.show()
    
    
    plt.savefig(output_PNG_filename)
    plt.close()
    print('\nThe figure has been successfully generated as: ', output_PNG_filename, '\n')
    

    
    executionTime = (datetime.now() - startTime)
    current_time = datetime.now().strftime("%Y/%m/%d at %H:%M:%S")
    print('\n\nThe script completed at: ' + str(current_time))
    print('Execution time: ' + str(executionTime), ' \n\n')

