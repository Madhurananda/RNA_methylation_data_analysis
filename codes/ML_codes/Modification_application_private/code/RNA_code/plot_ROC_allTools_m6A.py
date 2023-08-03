
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

from plot_ROC_allTools_func import *

'''
This script draws the ROC curves for the other existing tools. 
It reads the output from those tools. Tools are: 

1. Epinano
2. m6Anet
3. Tombo 
4. DRUMMER
5. xPore
6. CHEUI
7. Our Tool
'''



if __name__=='__main__' and '__file__' in globals():

    startTime = datetime.now()
    current_time = startTime.strftime("%Y/%m/%d at %H:%M:%S")
    print('\n\nThe script starting at: ' + str(current_time), ' \n\n' )

    ## This data token is consistent for one data type and all existing tools 
    base_dir = '/home/madhu/work/comp_tools/paper_6/download_tools'



    #####################################################################################
    '''
    This one is about IVT m6A data
    '''
    data_token = 'IVT_m6A'
    ref_file='/home/madhu/work/ref_transcriptome/IVT_seq/IVT_seq.fa'
    ref_content = read_file(ref_file)
    CHR_SELECTED = ['A1', 'A2']
    NCLD_SELECTED = ['A', 'T']

    drummer_dir = base_dir+'/DRUMMER/RESULTS'
    DRUMMER_wt = get_DRUMMER_files( drummer_dir+"/IVT_m6A_0/IVT_m6A_0__isoform/IVT_normalA_fast5_0-IVT_normalA_fast5_1/complete_analysis/*.txt", CHR_SELECTED)
    DRUMMER_ko = get_DRUMMER_files( drummer_dir+"/IVT_m6A_1/IVT_m6A_1__isoform/IVT_normalA_fast5_1-IVT_m6A_fast5/complete_analysis/*.txt", CHR_SELECTED)

    xPore_dir = base_dir+'/xPore/RESULTS'
    data_token_tmp = 'IVT_m6A_0'
    xPore_wt = glob.glob( xPore_dir+'/'+data_token_tmp+'out_diffmod/diffmod.table' )
    data_token_tmp = 'IVT_m6A_1'
    xPore_ko = glob.glob( xPore_dir+'/'+data_token_tmp+'out_diffmod/diffmod.table' )

    our_log_file = '/home/madhu/Logs/predict__train_p3_BOTH_pred_p5_m6A_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.240025.pt__A.log'    
    # our_log_file = '/home/madhu/Logs/predict__train_p3_BOTH_pred_p5_m6A_A_BilstmMean.layers3.hs1024.F.lr1500.b512.p1000.GPU.ep1.240025.pt__A.log'

    #####################################################################################

    #####################################################################################

    '''
    This one is about m6A data in paper 3 
    '''
    # data_token = 'm6A_p3__49_50'
    # ref_file='/home/madhu/work/ref_transcriptome/OTHERS/synthetic_construct/Best_Practices_dRNAseq_analysis/reference_fasta/cc.fasta'
    # ref_content = read_file(ref_file)

    # NCLD_SELECTED = ['A', 'T']
    # CHR_SELECTED = []
    # selected_chr = []
    # for a_part in (ref_content.split('\n')):
    #     if a_part.startswith('>'):
    #         # selected_chr.append( a_part.split('>')[1] )
    #         CHR_SELECTED.append( a_part.split('>')[1] )
    
    # print('CHR_SELECTED: ', CHR_SELECTED)

    # drummer_dir = base_dir+'/DRUMMER/RESULTS'
    # DRUMMER_wt = [file for file in glob.glob(drummer_dir+"/m6A_p3__49_51/m6A_p3__49_51__isoform/GSM3528749-GSM3528751/complete_analysis/*.txt")]
    # DRUMMER_ko = [file for file in glob.glob(drummer_dir+"/m6A_p3__49_50/m6A_p3__49_50__isoform/GSM3528749-GSM3528750/complete_analysis/*.txt")]
    # # DRUMMER_ko = [file for file in glob.glob(drummer_dir+"/m6A_p3__51_52/m6A_p3__51_52__isoform/GSM3528751-GSM3528752/complete_analysis/*.txt")]

    # xPore_dir = base_dir+'/xPore/RESULTS'
    # data_token_tmp = 'm6A_p3__49_51'
    # xPore_wt = glob.glob( xPore_dir+'/'+data_token_tmp+'out_diffmod/diffmod.table' )
    # data_token_tmp = 'm6A_p3__49_50'
    # xPore_ko = glob.glob( xPore_dir+'/'+data_token_tmp+'out_diffmod/diffmod.table' )

    # # our_log_file = '/home/madhu/Logs/predict__pred_p3_train_p5_m6A_ALL_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.270255.pt__A.log'
    # our_log_file = '/home/madhu/Logs/predict__pred_p3_train_p5_m6A_ALL_A_BilstmMean.layers3.hs1024.F.lr500.b64.p1000.GPU.ep1.310658.pt__A.log'

    #####################################################################################

    
    cheui_dir = base_dir+'/CHEUI/RESULTS'
    cheui_wt = glob.glob( cheui_dir+'/'+data_token+'/WT/WT_site_level_pred.txt' )
    cheui_ko = glob.glob( cheui_dir+'/'+data_token+'/KO/KO_site_level_pred.txt' )

    epinano_dir = base_dir+'/EpiNano/RESULTS'
    epinano_wt = glob.glob( epinano_dir+'/'+data_token+'/WT_'+data_token+'.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv' )
    epinano_ko = glob.glob( epinano_dir+'/'+data_token+'/KO_'+data_token+'.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv' )

    m6Anet_dir = base_dir+'/m6anet/RESULTS'
    m6Anet_wt = glob.glob( m6Anet_dir+'/'+data_token+'/results/WT/data.site_proba.csv' )
    m6Anet_ko = glob.glob( m6Anet_dir+'/'+data_token+'/results/KO/data.site_proba.csv' )

    tombo_dir = base_dir+'/Tombo_Output/RESULTS'
    tombo_wt = glob.glob( tombo_dir+'/'+data_token+'/'+data_token+'final_text__.fraction_modified_reads.m6A.plus.wig' )
    tombo_ko = glob.glob( tombo_dir+'/'+data_token+'/'+data_token+'final_text__.fraction_modified_reads.m6A.minus.wig' )





    ### DRUMMER
    df_ROC_DRUMMER = get_final_array('DRUMMER', DRUMMER_wt, DRUMMER_ko, ref_content, 'chr_ctrl', 'pos_ctrl', 'ref_ctrl', 'p_values_OR_adj', 0, 0, 1, CHR_SELECTED, NCLD_SELECTED)

    ### xPore 
    df_ROC_xPore = get_final_array('xPore', xPore_wt, xPore_ko, ref_content, 'id', 'position', 'kmer', 'pval_ko_vs_wt', 1, 5, 1, CHR_SELECTED, NCLD_SELECTED)

    ### CHEUI
    ## CHEUI states the position of 'CACTATAGC' as 30 and it starts at the position 31 (in my calculation) in the ref IVT seq: A1 
    ## They consider the 9-mer 
    df_ROC_CHEUI = get_final_array('CHEUI', cheui_wt, cheui_ko, ref_content, 'contig', 'position', 'site', 'probability', 5, 9, 0, CHR_SELECTED, NCLD_SELECTED)

    ### Epinano
    df_ROC_Epinano = get_final_array('Epinano', epinano_wt, epinano_ko, ref_content, 'Ref', 'position', '#Kmer', 'ProbM', 0, 5, 0, CHR_SELECTED, NCLD_SELECTED)

    ### m6Anet
    df_ROC_m6Anet = get_final_array('m6Anet', m6Anet_wt, m6Anet_ko, ref_content, 'transcript_id', 'transcript_position', 'kmer', 'probability_modified', 1, 5, 0, CHR_SELECTED, NCLD_SELECTED)

    ### Tombo
    df_ROC_Tombo = get_final_array('Tombo', tombo_wt, tombo_ko, ref_content, 0,0,0,0,0,0,0, CHR_SELECTED, NCLD_SELECTED)

    ### OUR DNN TOOL
    df_ROC_OUR = get_final_array('OUR', our_log_file, 0, ref_content, 0,0,0,0,0,0,0, CHR_SELECTED, NCLD_SELECTED)






    df_ROC = df_ROC_DRUMMER.join(df_ROC_xPore['xPore_probs__0']).join(df_ROC_xPore['xPore_probs__1']).join(df_ROC_CHEUI['CHEUI_probs__0']).join(df_ROC_CHEUI['CHEUI_probs__1']).join(df_ROC_Epinano['Epinano_probs__0']).join(df_ROC_Epinano['Epinano_probs__1']).join(df_ROC_m6Anet['m6Anet_probs__0']).join(df_ROC_m6Anet['m6Anet_probs__1']).join(df_ROC_Tombo['Tombo_probs__0']).join(df_ROC_Tombo['Tombo_probs__1']).join(df_ROC_OUR['OUR_probs__0']).join(df_ROC_OUR['OUR_probs__1'])


    print('df_ROC ::::::    ')
    print(df_ROC)

    save_CSV_file = '/home/madhu/work/comp_tools/paper_6/download_tools/ROC_outputs/roc_info_FINAL_'+data_token+'.csv'
    df_ROC.to_csv(save_CSV_file, index=False)

    




    # ## Load the CSV file: 
    df_ROC = pd.read_csv(save_CSV_file)

    # # print(df_ROC)



    y_true_DRUMMER, y_score_DRUMMER = get_true_scores(df_ROC, 'DRUMMER')
    y_true_xPore, y_score_xPore = get_true_scores(df_ROC, 'xPore')
    y_true_CHEUI, y_score_CHEUI = get_true_scores(df_ROC, 'CHEUI')
    y_true_epinano, y_score_epinano = get_true_scores(df_ROC, 'Epinano')
    y_true_m6Anet, y_score_m6Anet = get_true_scores(df_ROC, 'm6Anet')
    y_true_Tombo, y_score_Tombo = get_true_scores(df_ROC, 'Tombo')
    y_true_ours, y_score_ours = get_true_scores(df_ROC, 'OUR')





    # get_AUC(y_true_DRUMMER, y_score_DRUMMER)
    # get_AUC(y_true_xPore, y_score_xPore)
    # get_AUC(y_true_CHEUI, y_score_CHEUI)
    # get_AUC(y_true_epinano, y_score_epinano)
    # get_AUC(y_true_m6Anet, y_score_m6Anet)



    saveFig_title = base_dir + '/ROC_'+data_token+'.png'
    # values is a list containing: [y_true, y_score, tool]
    values = [y_true_ours, y_score_ours, 'Our_DNN_tool'], [y_true_epinano, y_score_epinano, 'Epinano'], [y_true_m6Anet, y_score_m6Anet, 'm6Anet'], [y_true_Tombo, y_score_Tombo, 'Tombo'], [y_true_DRUMMER, y_score_DRUMMER, 'DRUMMER'], [y_true_xPore, y_score_xPore, 'xPore'], [y_true_CHEUI, y_score_CHEUI, 'CHEUI']
    calc_AUC__draw_ROC(values, saveFig_title)



    executionTime = (datetime.now() - startTime)
    current_time = datetime.now().strftime("%Y/%m/%d at %H:%M:%S")
    print('\n\nThe script completed at: ' + str(current_time))
    print('Execution time: ' + str(executionTime), ' \n\n')


