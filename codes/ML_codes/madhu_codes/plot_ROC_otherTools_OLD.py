
import os,sys
# import time, datetime
sys.path.insert(0, '/home/madhu/work/codes')
from main import check_endWith, check_create_dir, check_env, read_file, write_file, count_lines, save_file

check_env('torch')

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

    fig, ax = plt.subplots(figsize=(20, 10))

    list_c = ['r', 'b', 'm']
    list_m = ['*', '.', 'o']

    # for a_value in values:
    for i in range(len(values)): 
        [y_true, y_score, tool] = values[i]
        fpr, tpr, thresholds = roc_curve(y_true, y_score)
        AUC_score = roc_auc_score(y_true, y_score)
        print('\n\nThe AUC is :', round(AUC_score, 4), ' for the TOOL: ', tool)

        plt.plot(fpr, tpr, list_c[i], list_m[i], linewidth=2, label='{}: ROC curve [Area Under Curve (AUC = {:.4f})]'.format(tool, AUC_score))
        plt.legend(loc='lower right', fontsize=20)
        plt.plot([0, 1], [0, 1], color='navy', linestyle='--')
        plt.xlim([-0.02, 1.02])
        plt.ylim([-0.02, 1.02])
        plt.xlabel('Specificity', fontsize=20)
        plt.ylabel('Sensitivity', fontsize=20)
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
    plt.xticks(np.linspace(0, 1, n_split), new_xTick_list, fontsize=15)
    yTick_list = []
    for n in np.linspace(0, 1, n_split):
        yTick_list.append(str(int(n*100))+'%')
    plt.yticks(np.linspace(0, 1, n_split), yTick_list, fontsize=15)
    plt.grid(color='y', linewidth=0.5)
    
    # plt.title('Mean ROC curve for TB Index Score', fontsize=35)
    plt.show()
    plt.savefig(saveFig_title)
    print('The ROC figure has been saved at: ', saveFig_title)
    plt.close('all')
    
    save_file(save_ROC_results, (fpr, tpr, thresholds, AUC_score))

    # return AUC_score



if __name__=='__main__' and '__file__' in globals():

    startTime = datetime.now()
    current_time = startTime.strftime("%Y/%m/%d at %H:%M:%S")
    print('\n\nThe script starting at: ' + str(current_time), ' \n\n' )

    ## This data token is consistent for one data type and all existing tools 
    base_dir = '/home/madhu/work/comp_tools/paper_6/download_tools'
    data_token = 'IVT_m6A'







    ### Epinano
    ## The output file to be read 
    epinano_dir = base_dir+'/EpiNano/RESULTS'
    out_files = glob.glob( epinano_dir+'/'+data_token+'/pred_DELTA*.csv' )
    ko_out_files = glob.glob( epinano_dir+'/'+data_token+'/KO_'+data_token+'.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv' )
    wt_out_files = glob.glob( epinano_dir+'/'+data_token+'/WT_'+data_token+'.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv' )
    
    epinano_df = read_outputs(out_files)
    ko_epinano_df = read_outputs(ko_out_files)
    wt_epinano_df = read_outputs(wt_out_files)

    # print(epinano_df)
    # print(epinano_df['ProbM'].values)
    # print(epinano_df['ProbU'].values)
    
    y_true_epinano = []
    y_score_epinano = []

    y_true_epinano.extend( len(ko_epinano_df)*[1] )
    y_score_epinano.extend( ko_epinano_df['ProbM'].values )

    y_true_epinano.extend( len(wt_epinano_df)*[0] )
    y_score_epinano.extend( wt_epinano_df['ProbM'].values )

    
    # save_ROC_results =  epinano_dir+'/'+data_token+'/Epinano_ROC_'+data_token+'results.pkl'
    # saveFig_title = epinano_dir+'/'+data_token+'/Epinano_ROC_'+data_token+'.png'
    # Epinano_AUC_score = draw_ROC(y_true, y_score, save_ROC_results, saveFig_title)
    # print('\n\nAUC score for Epinano: ', round(Epinano_AUC_score, 2))

    







    ### m6Anet
    m6Anet_dir = base_dir+'/m6anet/RESULTS'

    wt_out_files = glob.glob( m6Anet_dir+'/'+data_token+'/results/WT/data.site_proba.csv' )
    ko_out_files = glob.glob( m6Anet_dir+'/'+data_token+'/results/KO/data.site_proba.csv' )

    # out_files = glob.glob( m6Anet_dir+'/'+data_token+'/results/WT/data.indiv_proba.csv' )

    # m6Anet_df = read_outputs(out_files)
    # print(m6Anet_df)
    # print(m6Anet_df['probability_modified'].values)
    ko_m6Anet_df = read_outputs(ko_out_files)
    wt_m6Anet_df = read_outputs(wt_out_files)

    y_true_m6Anet = []
    y_score_m6Anet = []

    y_true_m6Anet.extend( len(ko_m6Anet_df)*[1] )
    y_score_m6Anet.extend( ko_m6Anet_df['probability_modified'].values )

    y_true_m6Anet.extend( len(wt_m6Anet_df)*[0] )
    y_score_m6Anet.extend( wt_m6Anet_df['probability_modified'].values )

    # save_ROC_results =  m6Anet_dir+'/'+data_token+'/results/m6Anet_ROC_'+data_token+'results.pkl'
    # saveFig_title = m6Anet_dir+'/'+data_token+'/results/m6Anet_ROC_'+data_token+'.png'
    # m6Anet_AUC_score = draw_ROC(y_true, y_score, save_ROC_results, saveFig_title)
    # print('\n\nAUC score for m6Anet: ', round(m6Anet_AUC_score, 2))








    
    ### Tombo

    tombo_dir = base_dir+'/Tombo_Output/RESULTS'
    # wt_out_files = glob.glob( tombo_dir+'/'+data_token+'/IVT_m6Afinal_text__.fraction_modified_reads.m6A.minus.wig' )
    # ko_out_files = glob.glob( tombo_dir+'/'+data_token+'/IVT_m6Afinal_text__.fraction_modified_reads.m6A.plus.wig' )

    ## wt_out_files = glob.glob( tombo_dir+'/'+data_token+'/IVT_m6Afinal_text__.fraction_modified_reads.m6A.plus.wig' )
    ## ko_out_files = glob.glob( tombo_dir+'/'+data_token+'/IVT_m6Afinal_text__.fraction_modified_reads.m6A.minus.wig' )

    wt_out_files = glob.glob( tombo_dir+'/'+data_token+'/IVT_m6Afinal_text__.fraction_modified_reads.m6A.plus.wig' )
    ko_out_files = glob.glob( tombo_dir+'/'+data_token+'/IVT_m6Afinal_text__.fraction_modified_reads.m6A.minus.wig' )

    ## wt_out_files = glob.glob( tombo_dir+'/'+data_token+'/IVT_m6Afinal_text__.dampened_fraction_modified_reads.m6A.minus.wig' )
    ## ko_out_files = glob.glob( tombo_dir+'/'+data_token+'/IVT_m6Afinal_text__.dampened_fraction_modified_reads.m6A.plus.wig' )

    # wt_out_files = glob.glob( tombo_dir+'/'+data_token+'/IVT_m6Afinal_text__.dampened_fraction_modified_reads.m6A.plus.wig' )
    # ko_out_files = glob.glob( tombo_dir+'/'+data_token+'/IVT_m6Afinal_text__.dampened_fraction_modified_reads.m6A.minus.wig' )
    

    y_true_Tombo = []
    y_score_Tombo = []
    y_true_Tombo, y_score_Tombo = read_Tombo_outs('wt', wt_out_files[0], y_true_Tombo, y_score_Tombo)
    y_true_Tombo, y_score_Tombo = read_Tombo_outs('ko', ko_out_files[0], y_true_Tombo, y_score_Tombo)

    # print('y_true: ', len(y_true))
    # print('y_score: ', len(y_score))

    # save_ROC_results =  tombo_dir+'/'+data_token+'/Tombo_ROC_'+data_token+'results.pkl'
    # saveFig_title = tombo_dir+'/'+data_token+'/Tombo_ROC_'+data_token+'.png'
    # Tombo_AUC_score = draw_ROC(y_true, y_score, save_ROC_results, saveFig_title)
    # print('\n\nAUC score for Tombo: ', round(Tombo_AUC_score, 2))






    ### DRUMMER
    drummer_dir = base_dir+'/DRUMMER/RESULTS'
    out_files = glob.glob( drummer_dir+'/'+data_token+'/'+ data_token+'__isoform/*/summary.txt' )

    # print('The file is: ', out_files[0])

    file_content = []
    for a_line in read_file(out_files[0]).split('\n'):
        if not len(a_line.split('\t')) == 1:
            file_content.append(a_line.split('\t'))


    # print(np.array(file_content))
    drummer_df = pd.DataFrame(data=np.array(file_content))
    drummer_df.columns = drummer_df.iloc[0]
    drummer_df = drummer_df.reset_index(drop=True)
    drummer_df = drummer_df[1:]
    
    # print(drummer_df['OR_padj'].values)
    # print(drummer_df['G_padj'].values)








    ### xPore 
    xPore_dir = base_dir+'/xPore/RESULTS'
    out_files = glob.glob( xPore_dir+'/'+data_token+'out_diffmod/diffmod.table' )

    # print('The file is: ', out_files[0])

    # xPore_df = read_outputs(out_files)
    # # print(xPore_df)
    # print(xPore_df['z_score_ko_vs_wt'].values)

    

    saveFig_title = base_dir + '/ROC_'+data_token+'.png'
    
    # values is a list containing: [y_true, y_score, tool]
    [y_true_epinano, y_score_epinano, 'Epinano'], [y_true_m6Anet, y_score_m6Anet, 'm6Anet'], [y_true_Tombo, y_score_Tombo, 'Tombo']
    calc_AUC__draw_ROC([y_true_epinano, y_score_epinano, 'Epinano'], [y_true_m6Anet, y_score_m6Anet, 'm6Anet'], [y_true_Tombo, y_score_Tombo, 'Tombo'], saveFig_title)


    
    executionTime = (datetime.now() - startTime)
    current_time = datetime.now().strftime("%Y/%m/%d at %H:%M:%S")
    print('\n\nThe script completed at: ' + str(current_time))
    print('Execution time: ' + str(executionTime), ' \n\n')


    