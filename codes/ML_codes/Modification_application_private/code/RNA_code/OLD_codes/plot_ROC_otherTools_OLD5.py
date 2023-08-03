
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

def get_values(result_file, ref_content, col_name):

    '''
    For m6A, only consider chr: A1 and A2. 
    '''
    print()
    print('The file is: ', result_file)
    file_content = []
    for a_line in read_file(result_file).split('\n'):
        if not len(a_line.split('\t')) == 1:
            file_content.append(a_line.split('\t'))

    # print(np.array(file_content))
    drummer_df = pd.DataFrame(data=np.array(file_content))
    drummer_df.columns = drummer_df.iloc[0]
    drummer_df = drummer_df.reset_index(drop=True)
    drummer_df = drummer_df[1:]

    # print(drummer_df)
    # print(len(list(drummer_df.columns)))


    selected_chr = drummer_df['chr_ctrl'].unique()
    # print(selected_chr)
    
    if len(selected_chr) > 1:
        print('*\n\nThere is something wrong with the selected chr. INVESTIGATE....')
    
    selected_chr_FINAL = selected_chr[0]

    if selected_chr_FINAL == 'A1' or selected_chr_FINAL == 'A2':

        selected_chr_seq = ref_content.split('>'+selected_chr_FINAL)[1].split('>')[0]
        print('The number of ALL nucleotides in the reference sequence: ', len(selected_chr_seq), ' for chromosome:', selected_chr_FINAL)



        list_chr = []
        list_seq_pos = []
        list_pVals = []
        chr_list = ['A', 'C', 'G', 'T']
        # pVal_dict_N_mthyl = {}
        for a_chr in chr_list:
            # if a_chr == 'C' or a_chr == 'G' or a_chr == 'T':
            if a_chr == 'C' or a_chr == 'G':
                continue

            pVal_list_N_mthyl = []

            seq_pos_list = natsorted([i for i, j in enumerate(selected_chr_seq) if j == a_chr])
            print('Chr: ', a_chr, ' appears ', len(seq_pos_list), ' times in the reference sequence.')

            drummer_df__chr = drummer_df[drummer_df['ref_ctrl'] == a_chr]
            print('Chr: ', a_chr, ' appears ', len(drummer_df__chr), ' times in the DRUMMER output.')
            # print( 'pos_list: ', seq_pos_list )
            # print( 'drummer_df__chr: ', pd.to_numeric(drummer_df__chr['pos_ctrl'].values) )
            pos_DIFF = [x for x in seq_pos_list if x not in pd.to_numeric(drummer_df__chr['pos_ctrl'].values) ]
            if len(pos_DIFF) > 0:
                print('The positions where ', a_chr, ' appears in the reference sequence, but not in the DRUMMER output: ', [(x-1) for x in pos_DIFF] )
            
            [pVal_list_N_mthyl.insert(x-1, 1) for x in pos_DIFF]
            # pVal_list_N_mthyl.extend( [1]*len(pos_DIFF) )
            
            for a_pos in seq_pos_list:
                for i in range(len(drummer_df__chr)):
                    
                    if a_pos == int(drummer_df__chr['pos_ctrl'].values[i]) :
                        # pVal_list_N_mthyl.append( np.float_(drummer_df__chr['padj'].values[i]) )
                        pVal_list_N_mthyl.append(-np.log( np.float_(drummer_df__chr[col_name].values[i]) if np.float_(drummer_df__chr[col_name].values[i])>10**(-300) else 10**(-300)) )

            list_chr.extend( [a_chr]*len(pVal_list_N_mthyl) )
            list_seq_pos.extend( seq_pos_list )
            list_pVals.extend( pVal_list_N_mthyl )

        final_array = np.transpose( np.vstack(( np.array([selected_chr_FINAL]*len(list_chr)), np.array(list_chr), np.array(list_seq_pos), np.array(list_pVals) )) )
        print('The shape of the final array: ', final_array.shape)
        print()
    else:
        final_array = 0

    return final_array

def file_content_array(out_files):
    for i in range(len(out_files)):
        array_out = get_values(out_files[i], ref_content, 'padj')
        if array_out != 0:
            if i == 0:
                array = array_out
            else:
                array = np.vstack(( array, array_out )) 
    return array


if __name__=='__main__' and '__file__' in globals():

    startTime = datetime.now()
    current_time = startTime.strftime("%Y/%m/%d at %H:%M:%S")
    print('\n\nThe script starting at: ' + str(current_time), ' \n\n' )

    ## This data token is consistent for one data type and all existing tools 
    base_dir = '/home/madhu/work/comp_tools/paper_6/download_tools'
    data_token = 'IVT_m6A'




    ## DRUMMER
    


    drummer_dir = base_dir+'/DRUMMER/RESULTS'
    y_true_DRUMMER = []
    y_score_DRUMMER = []

    ref_file='/home/madhu/work/ref_transcriptome/IVT_seq/IVT_seq.fa'
    ref_content = read_file(ref_file)

    # data_token_tmp = 'IVT_m6A_0'
    ## I need to sort this out later .... 
    # out_files = glob.glob( drummer_dir+'/'+data_token+'/'+ data_token+'__isoform/*/summary.txt' ) 
    # out_files = ['/home/madhu/work/comp_tools/paper_6/download_tools/DRUMMER/RESULTS/IVT_m6A_0/IVT_m6A_0__isoform/IVT_normalA_fast5_1-IVT_m6A_fast5_1/summary.txt']

    # out_files = ['/home/madhu/work/comp_tools/paper_6/download_tools/DRUMMER/RESULTS/IVT_m6A_0/IVT_m6A_0__isoform/IVT_normalA_fast5_0-IVT_normalA_fast5_1/complete_analysis/A1.complete.txt']
    out_files__0 = [file for file in glob.glob("/home/madhu/work/comp_tools/paper_6/download_tools/DRUMMER/RESULTS/IVT_m6A_0/IVT_m6A_0__isoform/IVT_normalA_fast5_0-IVT_normalA_fast5_1/complete_analysis/*.txt")]

    print()

    # print('The file is: ', out_files)

    neg_array = file_content_array(out_files__0)

    # for i in range(len(out_files__0)):
    #     array_out = get_values(out_files__0[i], ref_content, 'padj')
    #     if array_out != 0:
    #         if i == 0:
    #             neg_array = array_out
    #         else:
    #             neg_array = np.vstack(( neg_array, array_out )) 

    # neg_array = get_values(out_files[0], ref_content, 'padj')
    print('neg_array: ', neg_array)
    print('neg_array shape', neg_array.shape)

    # df_ROC_neg = pd.DataFrame(neg_array, columns=col_names)
    # print('df_ROC_neg: ', df_ROC_neg)
    
    # y_true_DRUMMER.extend( len(pVal_list_N_mthyl)*[0] )
    # y_score_DRUMMER.extend( pVal_list_N_mthyl )





    # out_files = ['/home/madhu/work/comp_tools/paper_6/download_tools/DRUMMER/RESULTS/IVT_m6A_1/IVT_m6A_1__isoform/IVT_normalA_fast5_1-IVT_m6A_fast5/complete_analysis/A1.complete.txt']
    out_files__1 = [file for file in glob.glob("/home/madhu/work/comp_tools/paper_6/download_tools/DRUMMER/RESULTS/IVT_m6A_1/IVT_m6A_1__isoform/IVT_normalA_fast5_1-IVT_m6A_fast5/complete_analysis/*.txt")]

    pos_array = file_content_array(out_files__1)

    # for i in range(len(out_files__1)):
    #     if i == 0:
    #         pos_array = get_values(out_files__1[i], ref_content, 'padj')
    #     else:
    #         pos_array = np.vstack(( pos_array, get_values(out_files__1[i], ref_content, 'padj') )) 

    # pos_array = get_values(out_files[0], ref_content, 'padj')

    ## check that the first three columns from the two array are the same .... 


    final_array_df = np.append(neg_array, pos_array[:, 3].reshape(len(pos_array[:, 1]), 1), 1)

    # print(final_array_df)
    col_names = ['chr', 'ref_seq_Nucl', 'position', 'DRUMMER_pVals__0', 'DRUMMER_pVals__1']
    df_ROC = pd.DataFrame(columns=col_names)

    df_ROC = pd.DataFrame(final_array_df, columns=col_names)
    print('df_ROC_pos: ', df_ROC)

    # y_true_DRUMMER.extend( len(pVal_list_mthyl)*[1] )
    # y_score_DRUMMER.extend( pVal_list_mthyl )

    # # y_true_DRUMMER.extend( len(drummer_df)*[0] )
    # # y_score_DRUMMER.extend( list( np.float_(drummer_df['frac_diff'].values)) )


    # # print(drummer_df['p_val'].values)
    # # print(y_true_DRUMMER)
    # # print( y_score_DRUMMER )
    # # AUC_score = roc_auc_score(y_true_DRUMMER, list(np.float_(y_score_DRUMMER)))
    # AUC_score = roc_auc_score(y_true_DRUMMER, y_score_DRUMMER)
    # print('\n\n**** AUC: ', AUC_score, '****\n\n')










    ### xPore 
    xPore_dir = base_dir+'/xPore/RESULTS'
    y_true_xPore = []
    y_score_xPore = []


    data_token_tmp = 'IVT_m6A_0'
    wt_out_files = glob.glob( xPore_dir+'/'+data_token_tmp+'out_diffmod/diffmod.table' )
    wt_xPore_df = read_outputs(wt_out_files)

    print(wt_xPore_df)
    # neg_array = get_values(out_files[0], ref_content, 'padj')




    # print(wt_xPore_df['z_score_ko_vs_wt'].values)
    y_true_xPore.extend( len(wt_xPore_df)*[0] )
    y_score_xPore.extend( wt_xPore_df['pval_ko_vs_wt'].values)
    # y_score_xPore.extend( list(np.float_(wt_xPore_df['z_score_ko_vs_wt'].values) ))


    data_token_tmp = 'IVT_m6A_1'
    ko_out_files = glob.glob( xPore_dir+'/'+data_token_tmp+'out_diffmod/diffmod.table' )
    ko_xPore_df = read_outputs(ko_out_files)
    # print(ko_xPore_df['z_score_ko_vs_wt'].values)
    y_true_xPore.extend( len(ko_xPore_df)*[1] )
    y_score_xPore.extend( ko_xPore_df['pval_ko_vs_wt'].values)

    # print(y_true_xPore)
    # print( y_score_xPore )
    # AUC_score = roc_auc_score(y_true_DRUMMER, list(np.float_(y_score_DRUMMER)))
    # AUC_score = roc_auc_score(y_true_xPore, y_score_xPore)
    # print('\n\n**** AUC: ', AUC_score, '****\n\n')












    # ### OUR DNN TOOL

    # y_true_ours = []
    # y_score_ours = []

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
    # y_score_ours_KO = []
    # y_score_ours_WT = []
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
    #                 y_score_ours_KO.append(int(conf_mat[1])/(int(conf_mat[0])+int(conf_mat[1])))
    #                 # y_true_ours.append(1)
    #                 y_score_ours_WT.append(int(conf_mat[3])/(int(conf_mat[2])+int(conf_mat[3])))

    # # print('The length of sites: ', len(y_score_ours))
    # N_BIN = 20
    # hist_KO,bins_KO=np.histogram(np.array(y_score_ours_KO),bins=np.linspace(0,1,N_BIN))
    # hist_WT,bins_WT=np.histogram(np.array(y_score_ours_WT),bins=np.linspace(0,1,N_BIN))
    # print('hist_KO: ', hist_KO)
    # print('hist_WT: ', hist_WT)
    # # print('sum(hist): ', sum(hist))
    # print('bins_KO: ', bins_KO)
    # print('bins_WT: ', bins_WT)

    # fig, ax = plt.subplots(figsize=(20, 20))
    # plt.clf()
    # # plt.hist(np.array(y_score_ours), bins=np.linspace(0,1,N_BIN), histtype='bar', label='y_score distribution')
    # plt.hist([ np.array(y_score_ours_KO) , np.array(y_score_ours_WT) ], bins=np.linspace(0,1,N_BIN), histtype='bar', label=['KO', 'WT'])
    # plt.legend(loc='best', fontsize=25)
    # plt.xticks([x for x in np.linspace(0,1,N_BIN)], [round(x, 2) for x in np.linspace(0,1,N_BIN)], fontsize=20)
    # plt.yticks(fontsize=20)
    # plt.grid(color='y', linestyle='-', linewidth=1)
    # plt.xlabel('Probablity Thresholds', fontsize=22)
    # plt.ylabel('Number of sites', fontsize=22)
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

