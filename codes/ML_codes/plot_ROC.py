
import os,sys
# import time, datetime
sys.path.insert(0, '/home/madhu/work/codes')
from main import check_endWith, check_create_dir, check_env, read_file, write_file, count_lines, save_file

import TS_finetune_pred_MP
from datetime import datetime

import warnings
# warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore")

from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib.pyplot as plt
import numpy as np

'''
This script draws the ROC curves.
At the moment, it reads the log output file and calculates the probs and true labels 
'''

if __name__=='__main__' and '__file__' in globals():

    check_env('torch')

    startTime = datetime.now()
    current_time = startTime.strftime("%Y/%m/%d at %H:%M:%S")
    print('\n\nThe script starting at: ' + str(current_time), ' \n\n' )

    # log_file_name = '/home/madhu/Logs/predict__train_p3_ALL_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.270011.pt__A.log'
    
    log_file_name = sys.argv[1]
    REF_GENOME = sys.argv[2]
    save_folder = sys.argv[3]
    check_endWith(save_folder)

    log_file_content = read_file(log_file_name)
    ref_gene_content = read_file(REF_GENOME)

    # pred_folder = sys.argv[9]
    # print(log_file_content.split('Pred: ')[0:2])
    # print(ref_gene_content)

    gene_list = []
    for a_part in ref_gene_content.split('\n'):
        if a_part.startswith('>'):
            gene_list.append( a_part.split('>')[1] )
            # print( a_part.split('>')[1] )
            # for a_gene_part in a_part.split('>'):
            #     print('a_gene_part: ', a_gene_part)
    y_true = []
    y_score = []
    for a_part in log_file_content.split('\n'):
        # if a_part.startswith('Pred: '+):
        #     a_part.split()
        for a_gene in gene_list:
            if a_part.startswith('Pred: '+a_gene):
                info_part = a_part.split(' ')
                if len(info_part) == 6:
                    conf_mat = info_part[2:len(info_part)]
                    
                    if (int(conf_mat[0])+int(conf_mat[1])) != 0 and (int(conf_mat[2])+int(conf_mat[3])) != 0:
                        y_true.append(0)
                        y_score.append(int(conf_mat[1])/(int(conf_mat[0])+int(conf_mat[1])))
                        # print('y_true: ', y_true)
                        # print('y_score: ', y_score)

                        y_true.append(1)
                        y_score.append(int(conf_mat[3])/(int(conf_mat[2])+int(conf_mat[3])))
                        # print('y_true: ', y_true)
                        # print('y_score: ', y_score)
                    
                # if a_part.startswith(a_gene):
                #     print( a_part.split(' ') )

        # for a_value in a_part.split(' '):
        #     if a_value != '':
        #         print('a_value: ', a_value)

    print('y_true: ', y_true)
    print('y_score: ', y_score)


    '''
    Draw the ROc curves .... 
    '''

    saveFig_title = save_folder+save_folder.split('/')[len(save_folder.split('/'))-2]+'__ROC_curve.png'

    fpr, tpr, thresholds = roc_curve(y_true, y_score)
    AUC_score = roc_auc_score(y_true, y_score)

    save_ROC_results = save_folder+save_folder.split('/')[len(save_folder.split('/'))-2]+'__ROC_results.pkl'
    save_file(save_ROC_results, (fpr, tpr, thresholds, AUC_score))

    print('XXXXXXX The AUC score is: ', round(AUC_score, 2))
    fig, ax = plt.subplots(figsize=(20, 10))
    plt.plot(fpr, tpr, 'r', linewidth=2, label='ROC curve [Area Under Curve (AUC = {:.4f})]'.format(AUC_score))
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

    
    # executionTime = (datetime.now() - startTime)
    # current_time = datetime.now().strftime("%Y/%m/%d at %H:%M:%S")
    # print('\n\nThe script completed at: ' + str(current_time))
    # print('Execution time: ' + str(executionTime), ' \n\n')



