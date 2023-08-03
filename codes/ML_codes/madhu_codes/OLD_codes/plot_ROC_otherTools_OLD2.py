
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

'''
This script draws the ROC curves for the other existing tools. 
It reads the output from those tools. Tools are: 

1. Epinano
2. m6Anet
3. Tombo 
4. DRUMMER
5. xPore

'''

def read_outputs(csv_files):
    if len(csv_files) != 1:
        print('******** There are multiple csv files. INVESTIGATE ... *********')
    else:
        csv_file = csv_files[0]
    
    # print('The final output file is: ', csv_file)
    df = pd.read_csv(csv_file)
    return df

def read_Tombo_outs(data_type, out_file, y_true, y_score):
    for a_line in read_file(out_file).split('\n'):
        if len(a_line) >=1:
            if a_line[0].isdigit():
                # print(a_line.split(' ')[-1])
                if data_type == 'wt': 
                    y_true.append(0)
                elif data_type == 'ko': 
                    y_true.append(1)
                y_score.append( float(a_line.split(' ')[-1]) )
    return y_true, y_score


# def draw_ROC(y_true, y_score, save_ROC_results, saveFig_title):
# values is a list containing: [y_true, y_score, tool]
def calc_AUC__draw_ROC(values, saveFig_title): 

    fig, ax = plt.subplots(figsize=(20, 20))

    list_c = ['r', 'b', 'm', 'y', 'g', 'c']
    list_m = ['*', 'x', 'o', '.', '*', 'x']
    list_s = ['-', '--', '-.', ':', ' ', '']
    # list_s = ['-', '--', '-.', ':', ' ', '', 'solid', 'dashed', 'dashdot', 'dotted']

    if len(values) != len(list_c) or len(values) != len(list_c):
        print('***** SOMETHING is WRONG, INVESTIGATE ... ******')
        sys.exit(0)
    
    # for a_value in values:
    for i in range(len(values)): 
        [y_true, y_score, tool] = values[i]
        fpr, tpr, thresholds = roc_curve(y_true, y_score)
        AUC_score = roc_auc_score(y_true, y_score)
        print('\n\nThe AUC is :', round(AUC_score, 4), ' for the TOOL: ', tool)
        plt.plot(fpr, tpr, list_c[i], marker=list_m[i], linestyle=list_s[i], linewidth=2, label='{}: (AUC = {:.4f})'.format(tool, AUC_score))

    plt.legend(loc='lower right', fontsize=28)
    plt.plot([0, 1], [0, 1], color='navy', linestyle='--')
    plt.xlim([-0.02, 1.02])
    plt.ylim([-0.02, 1.02])
    plt.xlabel('Specificity', fontsize=25)
    plt.ylabel('Sensitivity', fontsize=25)
        # plt.xlabel('False Positive Rate', fontsize=20)
        # plt.ylabel('True Positive Rate', fontsize=20)
    
    n_split = 11
    xTick_list = []
    for n in np.linspace(0, 1, n_split):
        xTick_list.append(str(int(n*100))+'%')
    # reversing the list
    new_xTick_list = []
    for i in xTick_list:
        new_xTick_list.insert(0, i)
    plt.xticks(np.linspace(0, 1, n_split), new_xTick_list, fontsize=22)
    yTick_list = []
    for n in np.linspace(0, 1, n_split):
        yTick_list.append(str(int(n*100))+'%')
    plt.yticks(np.linspace(0, 1, n_split), yTick_list, fontsize=22)
    plt.grid(color='y', linewidth=0.5)
    
    # plt.title('Mean ROC curve for TB Index Score', fontsize=35)
    plt.show()
    plt.savefig(saveFig_title)
    print('The ROC figure has been saved at: ', saveFig_title)
    plt.close('all')
    
    # save_file(save_ROC_results, (fpr, tpr, thresholds, AUC_score))

    # return AUC_score



if __name__=='__main__' and '__file__' in globals():

    startTime = datetime.now()
    current_time = startTime.strftime("%Y/%m/%d at %H:%M:%S")
    print('\n\nThe script starting at: ' + str(current_time), ' \n\n' )

    ## This data token is consistent for one data type and all existing tools 
    base_dir = '/home/madhu/work/comp_tools/paper_6/download_tools'
    data_token = 'IVT_m6A'







    # ### Epinano
    # ## The output file to be read 
    # epinano_dir = base_dir+'/EpiNano/RESULTS'
    # out_files = glob.glob( epinano_dir+'/'+data_token+'/pred_DELTA*.csv' )
    # ko_out_files = glob.glob( epinano_dir+'/'+data_token+'/KO_'+data_token+'.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv' )
    # wt_out_files = glob.glob( epinano_dir+'/'+data_token+'/WT_'+data_token+'.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv' )
    
    # epinano_df = read_outputs(out_files)
    # ko_epinano_df = read_outputs(ko_out_files)
    # wt_epinano_df = read_outputs(wt_out_files)

    # # print(epinano_df)
    # # print(epinano_df['ProbM'].values)
    # # print(epinano_df['ProbU'].values)
    
    # y_true_epinano = []
    # y_score_epinano = []

    # y_true_epinano.extend( len(ko_epinano_df)*[1] )
    # y_score_epinano.extend( ko_epinano_df['ProbM'].values )

    # y_true_epinano.extend( len(wt_epinano_df)*[0] )
    # y_score_epinano.extend( wt_epinano_df['ProbM'].values )

    







    # ### m6Anet
    # m6Anet_dir = base_dir+'/m6anet/RESULTS'

    # wt_out_files = glob.glob( m6Anet_dir+'/'+data_token+'/results/WT/data.site_proba.csv' )
    # ko_out_files = glob.glob( m6Anet_dir+'/'+data_token+'/results/KO/data.site_proba.csv' )

    # # out_files = glob.glob( m6Anet_dir+'/'+data_token+'/results/WT/data.indiv_proba.csv' )

    # # m6Anet_df = read_outputs(out_files)
    # # print(m6Anet_df)
    # # print(m6Anet_df['probability_modified'].values)
    # ko_m6Anet_df = read_outputs(ko_out_files)
    # wt_m6Anet_df = read_outputs(wt_out_files)

    # y_true_m6Anet = []
    # y_score_m6Anet = []

    # y_true_m6Anet.extend( len(ko_m6Anet_df)*[1] )
    # y_score_m6Anet.extend( ko_m6Anet_df['probability_modified'].values )

    # y_true_m6Anet.extend( len(wt_m6Anet_df)*[0] )
    # y_score_m6Anet.extend( wt_m6Anet_df['probability_modified'].values )










    
    # ### Tombo

    # tombo_dir = base_dir+'/Tombo_Output/RESULTS'
    # # wt_out_files = glob.glob( tombo_dir+'/'+data_token+'/IVT_m6Afinal_text__.fraction_modified_reads.m6A.minus.wig' )
    # # ko_out_files = glob.glob( tombo_dir+'/'+data_token+'/IVT_m6Afinal_text__.fraction_modified_reads.m6A.plus.wig' )

    # ## wt_out_files = glob.glob( tombo_dir+'/'+data_token+'/IVT_m6Afinal_text__.fraction_modified_reads.m6A.plus.wig' )
    # ## ko_out_files = glob.glob( tombo_dir+'/'+data_token+'/IVT_m6Afinal_text__.fraction_modified_reads.m6A.minus.wig' )

    # wt_out_files = glob.glob( tombo_dir+'/'+data_token+'/IVT_m6Afinal_text__.fraction_modified_reads.m6A.plus.wig' )
    # ko_out_files = glob.glob( tombo_dir+'/'+data_token+'/IVT_m6Afinal_text__.fraction_modified_reads.m6A.minus.wig' )

    # ## wt_out_files = glob.glob( tombo_dir+'/'+data_token+'/IVT_m6Afinal_text__.dampened_fraction_modified_reads.m6A.minus.wig' )
    # ## ko_out_files = glob.glob( tombo_dir+'/'+data_token+'/IVT_m6Afinal_text__.dampened_fraction_modified_reads.m6A.plus.wig' )

    # # wt_out_files = glob.glob( tombo_dir+'/'+data_token+'/IVT_m6Afinal_text__.dampened_fraction_modified_reads.m6A.plus.wig' )
    # # ko_out_files = glob.glob( tombo_dir+'/'+data_token+'/IVT_m6Afinal_text__.dampened_fraction_modified_reads.m6A.minus.wig' )
    

    # y_true_Tombo = []
    # y_score_Tombo = []
    # y_true_Tombo, y_score_Tombo = read_Tombo_outs('wt', wt_out_files[0], y_true_Tombo, y_score_Tombo)
    # y_true_Tombo, y_score_Tombo = read_Tombo_outs('ko', ko_out_files[0], y_true_Tombo, y_score_Tombo)

    # # print('y_true: ', len(y_true))
    # # print('y_score: ', len(y_score))








    ## DRUMMER
    drummer_dir = base_dir+'/DRUMMER/RESULTS'
    y_true_DRUMMER = []
    y_score_DRUMMER = []

    data_token_tmp = 'IVT_m6A_0'
    ## I need to sort this out later .... 
    # out_files = glob.glob( drummer_dir+'/'+data_token+'/'+ data_token+'__isoform/*/summary.txt' ) 
    out_files = ['/home/madhu/work/comp_tools/paper_6/download_tools/DRUMMER/RESULTS/IVT_m6A_0/IVT_m6A_0__isoform/IVT_normalA_fast5_1-IVT_m6A_fast5_1/summary.txt']

    out_files = ['/home/madhu/work/comp_tools/paper_6/download_tools/DRUMMER/RESULTS/IVT_m6A_0/IVT_m6A_0__isoform/IVT_normalA_fast5_0-IVT_normalA_fast5_1/complete_analysis/A1.complete.txt']



    # print('The file is: ', out_files)

    file_content = []
    for a_line in read_file(out_files[0]).split('\n'):
        if not len(a_line.split('\t')) == 1:
            file_content.append(a_line.split('\t'))

    # print(np.array(file_content))
    drummer_df = pd.DataFrame(data=np.array(file_content))
    drummer_df.columns = drummer_df.iloc[0]
    drummer_df = drummer_df.reset_index(drop=True)
    drummer_df = drummer_df[1:]

    # print(drummer_df)
    # print(len(list(drummer_df.columns)))
    '''
    ['chr_ctrl', 'pos_ctrl', 'ref_ctrl', 'depth_ctrl', 'A_ctrl', 'C_ctrl', 'G_ctrl', 'T_ctrl', 'N_ctrl', 'ref_fraction_ctrl', 'chr_treat', 'ref_treat', 'depth_treat', 'A_treat', 'C_treat', 'G_treat', 'T_treat', 'N_treat', 'ref_fraction_treat', 'ratio_treat', 'frac_diff', 'ratio_ctrl', 'log2_(OR)', 'odds_ratio', 'p_values_OR', 'nearest_ac', 'nearest_ac_motif', 'five_bp_motif', 'eleven_bp_motif', 'G_test', 'p_val', 'padj', 'p_values_OR_adj', 'is_SNP', 'accumulation', 'depletion', 'homopolymer']
    '''
    for a_c in list(drummer_df.columns):
        print('\n\n*** ', a_c, ' : ')
        print(drummer_df[a_c].values)

    # print(drummer_df['p_val'].values)
    # print(drummer_df['OR_padj'].values)
    # print(drummer_df['G_padj'].values)
    y_true_DRUMMER.extend( len(drummer_df)*[0] )
    y_score_DRUMMER.extend( list(np.float_(drummer_df['padj'].values) ))

    # out_files = ['/home/madhu/work/comp_tools/paper_6/download_tools/DRUMMER/RESULTS/IVT_m6A_0/IVT_m6A_0__isoform/IVT_normalA_fast5_0-IVT_normalA_fast5_1/complete_analysis/A2.complete.txt']
    # file_content = []
    # for a_line in read_file(out_files[0]).split('\n'):
    #     if not len(a_line.split('\t')) == 1:
    #         file_content.append(a_line.split('\t'))
    # drummer_df = pd.DataFrame(data=np.array(file_content))
    # drummer_df.columns = drummer_df.iloc[0]
    # drummer_df = drummer_df.reset_index(drop=True)
    # drummer_df = drummer_df[1:]
    # print(drummer_df['p_val'].values)
    # # print(drummer_df['OR_padj'].values)
    # # print(drummer_df['G_padj'].values)
    # y_true_DRUMMER.extend( len(drummer_df)*[1] )
    # y_score_DRUMMER.extend( list(np.float_(drummer_df['p_val'].values) ))





    # data_token_tmp = 'IVT_m6A_1'
    # # out_files = glob.glob( drummer_dir+'/'+data_token+'/'+ data_token+'__isoform/*/summary.txt' ) 
    # # out_files = ['/home/madhu/work/comp_tools/paper_6/download_tools/DRUMMER/RESULTS/IVT_m6A_1/IVT_m6A_1__isoform/IVT_normalA_fast5_1-IVT_m6A_fast5/summary.txt']
    # out_files = ['/home/madhu/work/comp_tools/paper_6/download_tools/DRUMMER/RESULTS/IVT_m6A_1/IVT_m6A_1__isoform/IVT_normalA_fast5_1-IVT_m6A_fast5/complete_analysis/A1.complete.txt']

    # file_content = []
    # for a_line in read_file(out_files[0]).split('\n'):
    #     if not len(a_line.split('\t')) == 1:
    #         file_content.append(a_line.split('\t'))

    # # print(np.array(file_content))
    # drummer_df = pd.DataFrame(data=np.array(file_content))
    # drummer_df.columns = drummer_df.iloc[0]
    # drummer_df = drummer_df.reset_index(drop=True)
    # drummer_df = drummer_df[1:]
    
    # # print(drummer_df['OR_padj'].values)
    # # print(drummer_df['G_padj'].values)

    # y_true_DRUMMER.extend( len(drummer_df)*[1] )
    # y_score_DRUMMER.extend( list( np.float_(drummer_df['padj'].values)) )

 

    # # out_files = ['/home/madhu/work/comp_tools/paper_6/download_tools/DRUMMER/RESULTS/IVT_m6A_1/IVT_m6A_1__isoform/IVT_normalA_fast5_1-IVT_m6A_fast5/complete_analysis/A2.complete.txt']
    # # file_content = []
    # # for a_line in read_file(out_files[0]).split('\n'):
    # #     if not len(a_line.split('\t')) == 1:
    # #         file_content.append(a_line.split('\t'))
    # # drummer_df = pd.DataFrame(data=np.array(file_content))
    # # drummer_df.columns = drummer_df.iloc[0]
    # # drummer_df = drummer_df.reset_index(drop=True)
    # # drummer_df = drummer_df[1:]
    # # y_true_DRUMMER.extend( len(drummer_df)*[0] )
    # # y_score_DRUMMER.extend( list( np.float_(drummer_df['frac_diff'].values)) )


    # # print(drummer_df['p_val'].values)
    # # print(y_true_DRUMMER)
    # # print( y_score_DRUMMER )
    # # AUC_score = roc_auc_score(y_true_DRUMMER, list(np.float_(y_score_DRUMMER)))
    # AUC_score = roc_auc_score(y_true_DRUMMER, y_score_DRUMMER)
    # print('\n\n**** AUC: ', AUC_score, '****\n\n')










    # ### xPore 
    # xPore_dir = base_dir+'/xPore/RESULTS'
    # y_true_xPore = []
    # y_score_xPore = []


    # data_token_tmp = 'IVT_m6A_0'
    # wt_out_files = glob.glob( xPore_dir+'/'+data_token_tmp+'out_diffmod/diffmod.table' )
    # wt_xPore_df = read_outputs(wt_out_files)
    # # print(wt_xPore_df['z_score_ko_vs_wt'].values)
    # y_true_xPore.extend( len(wt_xPore_df)*[0] )
    # y_score_xPore.extend( wt_xPore_df['pval_ko_vs_wt'].values)
    # # y_score_xPore.extend( list(np.float_(wt_xPore_df['z_score_ko_vs_wt'].values) ))


    # data_token_tmp = 'IVT_m6A_1'
    # ko_out_files = glob.glob( xPore_dir+'/'+data_token_tmp+'out_diffmod/diffmod.table' )
    # ko_xPore_df = read_outputs(ko_out_files)
    # # print(ko_xPore_df['z_score_ko_vs_wt'].values)
    # y_true_xPore.extend( len(ko_xPore_df)*[1] )
    # y_score_xPore.extend( ko_xPore_df['pval_ko_vs_wt'].values)

    # print(y_true_xPore)
    # print( y_score_xPore )
    # # AUC_score = roc_auc_score(y_true_DRUMMER, list(np.float_(y_score_DRUMMER)))
    # AUC_score = roc_auc_score(y_true_xPore, y_score_xPore)
    # print('\n\n**** AUC: ', AUC_score, '****\n\n')


    









    ### OUR DNN TOOL

    y_true_ours = []
    y_score_ours = []

    # ## This is using the log file ..................................
    # log_file_name = '/home/madhu/Logs/predict__train_p3_ALL_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.190363.pt__A.log'
    # # log_file_name = '/home/madhu/Logs/predict__train_p3_ALL_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.270011.pt__A.log'
    # # log_file_name = '/home/madhu/Logs/predict__train_p3_BOTH_pred_p6_m6A_A_BilstmMean.layers3.hs1024.F.lr1500.b512.p1000.GPU.ep1.160151.pt__A_NEW.log'
    # REF_GENOME = '/mnt/labshare/share/reference_genome/EpiNano_Reference_sequences/cc.fasta'


    # # log_file_name = '/home/madhu/Logs/predict__pred_p6_train_p5_m6A_ALL_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.330287.pt__A_NEW.log'
    # # REF_GENOME = '/home/madhu/work/ref_transcriptome/mouse/mrna.fa'
    
    # log_file_content = read_file(log_file_name)
    # ref_gene_content = read_file(REF_GENOME)

    # gene_list = []
    # for a_part in tqdm(ref_gene_content.split('\n')):
    #     if a_part.startswith('>'):
    #         gene_list.append( a_part.split('>')[1] )
    # for a_part in tqdm(log_file_content.split('\n')):
    #     for a_gene in gene_list:
    #         if a_part.startswith('Pred: '+a_gene):
    #             info_part = a_part.split(' ')
    #             if len(info_part) == 6:
    #                 conf_mat = info_part[2:len(info_part)]
                    
    #                 if (int(conf_mat[0])+int(conf_mat[1])) != 0 and (int(conf_mat[2])+int(conf_mat[3])) != 0:
    #                     y_true_ours.append(0)
    #                     y_score_ours.append(int(conf_mat[1])/(int(conf_mat[0])+int(conf_mat[1])))
    #                     y_true_ours.append(1)
    #                     y_score_ours.append(int(conf_mat[3])/(int(conf_mat[2])+int(conf_mat[3])))

    # AUC_score = roc_auc_score(y_true_ours, y_score_ours)
    # print('\n\n**** AUC: ', AUC_score, '****\n\n')



    # ## This is using the saved pickle file .................................. (for real data: mouse)

    # pred_results = load_file('/home/madhu/work/class_outputs/p5_decoding_epitranscriptional_native_RNA_seq/model_PRED_output__pred_p6_train_p5_m6A_ALL_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.330287.pt__A/pred_results.pkl')

    # for mapc in pred_results:
    #     # print("Total sites for <{}>: +{} -{}\n".format( mapc, (len(pred_results[mapc]['+']) if '+' in pred_results[mapc] else 0), (len(pred_results[mapc]['-']) if '-' in pred_results[mapc] else 0) ) );
    #     for mapstn in pred_results[mapc]:
    #         for ref_pos in pred_results[mapc][mapstn]:
                
    #             t_pred = pred_results[mapc][mapstn][ ref_pos ]
    #             # if ref_pos == 2161:
    #             #     print('t_pred: ', t_pred)
                
    #             # print("Pred: {}:{}:{} {} {} {} {}".format( mapc, mapstn, ref_pos, t_pred[0][0], t_pred[0][1], t_pred[1][0], t_pred[1][1] ) )
    #             conf_mat = [t_pred[0][0], t_pred[0][1], t_pred[1][0], t_pred[1][1]]
    #             if (int(conf_mat[0])+int(conf_mat[1])) != 0 and (int(conf_mat[2])+int(conf_mat[3])) != 0:
    #                 # y_true_ours.append(0)
    #                 # y_score_ours.append(int(conf_mat[1])/(int(conf_mat[0])+int(conf_mat[1])))
    #                 y_true_ours.append(1)
    #                 y_score_ours.append(int(conf_mat[3])/(int(conf_mat[2])+int(conf_mat[3])))

    # print('The length of sites: ', len(y_score_ours))
    # N_BIN = 10
    # hist,bins=np.histogram(np.array(y_score_ours),bins=np.linspace(0,1,N_BIN))
    # print('hist: ', hist)
    # # print('sum(hist): ', sum(hist))
    # print('bins: ', bins)

    # fig, ax = plt.subplots(figsize=(20, 20))
    # plt.clf()
    # plt.hist(np.array(y_score_ours), bins=np.linspace(0,1,N_BIN), histtype='bar', label='y_score distribution')
    # plt.legend(loc='best', fontsize=25)
    # plt.xticks([x for x in np.linspace(0,1,N_BIN)], [round(x, 2) for x in np.linspace(0,1,N_BIN)], fontsize=20)
    # plt.yticks(fontsize=20)
    # plt.grid(color='y', linestyle='-', linewidth=1)
    # plt.xlabel('Probablity Thresholds', fontsize=22)
    # plt.ylabel('Number of y-scores', fontsize=22)
    # plt.show()
    
    # output_PNG_filename = '/home/madhu/work/class_outputs/p5_decoding_epitranscriptional_native_RNA_seq/model_PRED_output__pred_p6_train_p5_m6A_ALL_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.330287.pt__A/yscore_dist.png'
    # plt.savefig(output_PNG_filename)
    # plt.close()
    # print('\nThe figure has been successfully generated as: ', output_PNG_filename, '\n')
    







    # saveFig_title = base_dir + '/ROC_'+data_token+'.png'
    # # values is a list containing: [y_true, y_score, tool]
    # values = [y_true_ours, y_score_ours, 'Our_DNN_tool'], [y_true_epinano, y_score_epinano, 'Epinano'], [y_true_m6Anet, y_score_m6Anet, 'm6Anet'], [y_true_Tombo, y_score_Tombo, 'Tombo'], [y_true_DRUMMER, y_score_DRUMMER, 'DRUMMER'], [y_true_xPore, y_score_xPore, 'xPore']
    # calc_AUC__draw_ROC(values, saveFig_title)


    
    executionTime = (datetime.now() - startTime)
    current_time = datetime.now().strftime("%Y/%m/%d at %H:%M:%S")
    print('\n\nThe script completed at: ' + str(current_time))
    print('Execution time: ' + str(executionTime), ' \n\n')

