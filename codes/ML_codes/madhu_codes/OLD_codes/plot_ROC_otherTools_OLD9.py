
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



# def read_Tombo_outs_OLD(data_type, out_file, y_true, y_score):
#     for a_line in read_file(out_file).split('\n'):
#         if len(a_line) >=1:
#             if a_line[0].isdigit():
#                 print(a_line.split(' ')[0], a_line.split(' ')[1])
#                 # print(a_line.split(' ')[-1])
#                 if data_type == 'wt': 
#                     y_true.append(0)
#                 elif data_type == 'ko': 
#                     y_true.append(1)
#                 y_score.append( float(a_line.split(' ')[-1]) )
#     return y_true, y_score



def read_Tombo_outs(data_type, out_file, y_true, y_score):
    for a_line in read_file(out_file).split('\n'):
        if len(a_line) >=1:
            if a_line[0].isdigit():
                print(a_line.split(' ')[0], a_line.split(' ')[1])
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



def file_content_array(out_files):
    for i in range(len(out_files)):
        array_out = get_values(out_files[i], ref_content, 'padj')
        if array_out != 0:
            if i == 0:
                array = array_out
            else:
                array = np.vstack(( array, array_out )) 
    return array


def pVal_prob(p_val):
    # -np.log( np.float_(drummer_df__chr[col_name].values[i]) if np.float_(drummer_df__chr[col_name].values[i])>10**(-300) else 10**(-300))
    prob = -np.log( np.float_(p_val) if np.float_(p_val)>10**(-300) else 10**(-300))
    return prob


def get_values(result_file, ref_content, col_name):

    '''
    For m6A, only consider chr: A1 and A2. 
    '''
    print()
    print('The file is: ', result_file)


    ########## DRUMMER 
    file_content = []
    for a_line in read_file(result_file).split('\n'):
        if not len(a_line.split('\t')) == 1:
            file_content.append(a_line.split('\t'))

    drummer_df = pd.DataFrame(data=np.array(file_content))
    drummer_df.columns = drummer_df.iloc[0]
    drummer_df = drummer_df.reset_index(drop=True)
    drummer_df = drummer_df[1:]
    selected_chr = drummer_df['chr_ctrl'].unique()
    
    if len(selected_chr) > 1:
        print('*\n\nThere is something wrong with the selected chr. INVESTIGATE....')
    
    selected_chr_FINAL = selected_chr[0]

    if selected_chr_FINAL == 'A1' or selected_chr_FINAL == 'A2':

        selected_chr_seq = ref_content.split('>'+selected_chr_FINAL)[1].split('>')[0]
        print('The number of ALL nucleotides in the reference sequence: ', len(selected_chr_seq), ' for chromosome:', selected_chr_FINAL)

        list_chr = []
        list_seq_pos = []
        list_pVals = []
        ncld_list = ['A', 'C', 'G', 'T']
        for a_ncld in ncld_list:
            # if a_chr == 'C' or a_chr == 'G' or a_chr == 'T':
            if a_ncld == 'C' or a_ncld == 'G':
                continue

            pVal_list_N_mthyl = []

            seq_pos_list = natsorted([i for i, j in enumerate(selected_chr_seq) if j == a_ncld])
            print('Nucleotide: ', a_ncld, ' appears ', len(seq_pos_list), ' times in the reference sequence.')

            drummer_df__chr = drummer_df[drummer_df['ref_ctrl'] == a_ncld]
            print('Nucleotide: ', a_ncld, ' appears ', len(drummer_df__chr), ' times in the DRUMMER output.')
            # print( 'pos_list: ', seq_pos_list )
            # print( 'drummer_df__chr: ', pd.to_numeric(drummer_df__chr['pos_ctrl'].values) )
            pos_DIFF = [x for x in seq_pos_list if x not in pd.to_numeric(drummer_df__chr['pos_ctrl'].values) ]
            if len(pos_DIFF) > 0:
                print('The positions where ', a_ncld, ' appears in the reference sequence, but not in the DRUMMER output: ', [(x-1) for x in pos_DIFF] )
            
            [pVal_list_N_mthyl.insert(x-1, pVal_prob(1)) for x in pos_DIFF]
            # pVal_list_N_mthyl.extend( [1]*len(pos_DIFF) )
            
            for a_pos in seq_pos_list:
                for i in range(len(drummer_df__chr)):
                    
                    if a_pos == int(drummer_df__chr['pos_ctrl'].values[i]) :
                        # pVal_list_N_mthyl.append( np.float_(drummer_df__chr['padj'].values[i]) )
                        # pVal_list_N_mthyl.append( -np.log( np.float_(drummer_df__chr[col_name].values[i]) if np.float_(drummer_df__chr[col_name].values[i])>10**(-300) else 10**(-300)) )
                        pVal_list_N_mthyl.append( pVal_prob( drummer_df__chr[col_name].values[i] ) )

            list_chr.extend( [a_ncld]*len(pVal_list_N_mthyl) )
            list_seq_pos.extend( seq_pos_list )
            list_pVals.extend( pVal_list_N_mthyl )

        final_array = np.transpose( np.vstack(( np.array([selected_chr_FINAL]*len(list_chr)), np.array(list_chr), np.array(list_seq_pos), np.array(list_pVals) )) )
        print('The shape of the final array: ', final_array.shape)
        
    else:
        print('For m6A: only A1 and A2 are considered.')
        final_array = 0
    
    print()
    return final_array



def get_values_xPore(result_file, ref_content, col_name):

    print()
    print('The file is: ', result_file)

    ########### xPore 
    xPore_df = read_outputs(result_file)
    xPore_df = xPore_df.astype({'position':'int'})

    selected_chr_temp = xPore_df['id'].unique()
    # print(selected_chr)

    print('For m6A: only A1 and A2 are considered. These chromosomes are overlooked: ', [x for x in selected_chr_temp if x != 'A1' and x != 'A2'])

    selected_chr = []
    selected_chr.extend([x for x in selected_chr_temp if x == 'A1' or x == 'A2'])

    final_array = 0
    for selected_chr_FINAL in natsorted(selected_chr):
        print()
        
        selected_chr_seq = ref_content.split('>'+selected_chr_FINAL)[1].split('>')[0]
        print('The number of ALL nucleotides in the reference sequence: ', len(selected_chr_seq), ' for chromosome:', selected_chr_FINAL)
        
        
        xPore_df__selected_temp = xPore_df[xPore_df['id'] == selected_chr_FINAL]
        xPore_df__selected = xPore_df__selected_temp.reset_index()
        xPore_df__selected = xPore_df__selected.drop('index', axis=1)
        # print(xPore_df__selected)
        # print('Duplicated:     ', xPore_df__selected[xPore_df__selected.duplicated('position')] )

        # print(xPore_df[xPore_df['position'] == 21])
        # print(xPore_df__selected[xPore_df__selected['position'] == 21])




        list_chr = []
        list_seq_pos = []
        list_pVals = []
        ncld_list = ['A', 'C', 'G', 'T']
        # pVal_dict_N_mthyl = {}
        for a_ncld in ncld_list:
            # if a_chr == 'C' or a_chr == 'G' or a_chr == 'T':
            if a_ncld == 'C' or a_ncld == 'G':
                continue

            pVal_list_N_mthyl = []

            seq_pos_list = natsorted([i for i, j in enumerate(selected_chr_seq) if j == a_ncld])
            print('Nucleotide: ', a_ncld, ' appears ', len(seq_pos_list), ' times in the reference sequence.')

            
            # xPore_df__chr = xPore_df[xPore_df['kmer'] == a_chr]
            ind_l = []
            for i in range(len(xPore_df__selected)):
                if xPore_df__selected._get_value(i, 'kmer')[2] == a_ncld:
                    ind_l.append(i)
                    # print(xPore_df__selected._get_value(i, 'kmer'))
            xPore_df__chr = xPore_df__selected.loc[ind_l, :]

            print('Nucleotide: ', a_ncld, ' appears ', len(xPore_df__chr), ' times in the xPore output 5-MERS.')
            # print(xPore_df__chr[xPore_df__chr['position'] == 21])
            # # print(xPore_df__chr[xPore_df__chr.duplicated('position')] )
            # print('ALL dupl', len(xPore_df__chr[xPore_df__chr.duplicated()]) )
            # print('ALL dupl pos: ', len(xPore_df__chr[xPore_df__chr.duplicated('position')]) )

            xPore_df__chr = xPore_df__chr.drop_duplicates()

            # xPore_df__chr['position'].values
            new_pos_values = [x+1 for x in xPore_df__chr['position'].values]
            xPore_df__chr['position'] = new_pos_values

            xPore_df__chr = xPore_df__chr.sort_values(by=['position'])

            # print(xPore_df__chr[xPore_df__chr['position'] == 21])
            # print('ALL dupl', len(xPore_df__chr[xPore_df__chr.duplicated()]) )
            # print('ALL dupl pos: ', len(xPore_df__chr[xPore_df__chr.duplicated('position')]) )

            # print( 'pos_list: ', seq_pos_list )
            # print( 'xPore_df__chr: ', pd.to_numeric(xPore_df__chr['position'].values) )
            pos_DIFF = [x for x in seq_pos_list if x not in pd.to_numeric(xPore_df__chr['position'].values) ]
            # pos_DIFF = [x for x in pd.to_numeric(xPore_df__chr['position'].values) if x not in seq_pos_list ]
            if len(pos_DIFF) > 0:
                print('The positions where ', a_ncld, ' appears in the reference sequence, but not in the xPore output: ', [x for x in pos_DIFF] )
                print('Total mis-match: ', len(pos_DIFF))
            
            [pVal_list_N_mthyl.insert(x-1, pVal_prob(1)) for x in pos_DIFF]
            # pVal_list_N_mthyl.extend( [1]*len(pos_DIFF) )
            
            for a_pos in seq_pos_list:
                for i in range(len(xPore_df__chr)):
                    
                    if a_pos == int(xPore_df__chr['position'].values[i]) :
                        # pVal_list_N_mthyl.append( np.float_(drummer_df__chr['padj'].values[i]) )
                        # pVal_list_N_mthyl.append(-np.log( np.float_(xPore_df__chr[col_name].values[i]) if np.float_(xPore_df__chr[col_name].values[i])>10**(-300) else 10**(-300)) )
                        pVal_list_N_mthyl.append( pVal_prob( xPore_df__chr[col_name].values[i] ) )
                    # else:
                    #     print('Here, the position does not match.', a_pos)

            list_chr.extend( [a_ncld]*len(pVal_list_N_mthyl) )
            list_seq_pos.extend( seq_pos_list )
            list_pVals.extend( pVal_list_N_mthyl )
            # print(len(list_chr))
            # print(len(list_seq_pos))
            # print(len(list_pVals) )

        array_out = np.transpose( np.vstack(( np.array([selected_chr_FINAL]*len(list_chr)), np.array(list_chr), np.array(list_seq_pos), np.array(list_pVals) )) )
        print('The shape of the final array: ', array_out.shape, ' for the chr: ', selected_chr_FINAL)

        if final_array == 0:
            final_array = array_out
        else:
            final_array = np.vstack(( final_array, array_out )) 
        
        # print('The shape of the final array at the END: ', final_array.shape, ' at chr: ', a_ncld)

    print()
    return final_array


# def get_values_CHEUI(result_file, ref_content, col_name):

#     print()
#     print('The file is: ', result_file)

#     cheui_df = pd.read_csv(result_file[0], delimiter='\t')
#     cheui_df = cheui_df.astype({'position':'int'})

#     selected_chr_temp = cheui_df['contig'].unique()

#     print('For m6A: only A1 and A2 are considered. These chromosomes are overlooked: ', [x for x in selected_chr_temp if x != 'A1' and x != 'A2'])

#     selected_chr = []
#     selected_chr.extend([x for x in selected_chr_temp if x == 'A1' or x == 'A2'])

#     final_array = 0
#     for selected_chr_FINAL in natsorted(selected_chr):
#         print()
        
#         selected_chr_seq = ref_content.split('>'+selected_chr_FINAL)[1].split('>')[0]
#         print('The number of ALL nucleotides in the reference sequence: ', len(selected_chr_seq), ' for chromosome:', selected_chr_FINAL)
        
        
#         CHEUI_df__selected_temp = cheui_df[cheui_df['contig'] == selected_chr_FINAL]
#         CHEUI_df__selected = CHEUI_df__selected_temp.reset_index()
#         CHEUI_df__selected = CHEUI_df__selected.drop('index', axis=1)


#         list_chr = []
#         list_seq_pos = []
#         list_pVals = []
#         ncld_list = ['A', 'C', 'G', 'T']
#         for a_ncld in ncld_list:
#             # if a_chr == 'C' or a_chr == 'G' or a_chr == 'T':
#             if a_ncld == 'C' or a_ncld == 'G':
#                 continue

#             pVal_list_N_mthyl = []

#             seq_pos_list = natsorted([i for i, j in enumerate(selected_chr_seq) if j == a_ncld])
#             print('Nucleotide: ', a_ncld, ' appears ', len(seq_pos_list), ' times in the reference sequence.')

            
#             ind_l = []
#             for i in range(len(CHEUI_df__selected)):
#                 if CHEUI_df__selected._get_value(i, 'site')[4] == a_ncld:
#                     ind_l.append(i)
#             CHEUI_df__chr = CHEUI_df__selected.loc[ind_l, :]

#             print('Nucleotide: ', a_ncld, ' appears ', len(CHEUI_df__chr), ' times in the CHEUI output 9-MERS.')

#             CHEUI_df__chr = CHEUI_df__chr.drop_duplicates()
#             new_pos_values = [x+5 for x in CHEUI_df__chr['position'].values]
#             CHEUI_df__chr['position'] = new_pos_values
#             CHEUI_df__chr = CHEUI_df__chr.sort_values(by=['position'])

#             # print( 'pos_list: ', seq_pos_list )
#             # print( 'xPore_df__chr: ', pd.to_numeric(CHEUI_df__chr['position'].values) )
#             pos_DIFF = [x for x in seq_pos_list if x not in pd.to_numeric(CHEUI_df__chr['position'].values) ]

#             if len(pos_DIFF) > 0:
#                 print('The positions where ', a_ncld, ' appears in the reference sequence, but not in the CHEUI output: ', [x for x in pos_DIFF] )
#                 print('Total mis-match: ', len(pos_DIFF))
            
#             [pVal_list_N_mthyl.insert(x-1, pVal_prob(1)) for x in pos_DIFF]
            
#             for a_pos in seq_pos_list:
#                 for i in range(len(CHEUI_df__chr)):
                    
#                     if a_pos == int(CHEUI_df__chr['position'].values[i]) :
#                         pVal_list_N_mthyl.append(  CHEUI_df__chr[col_name].values[i] ) 

#             list_chr.extend( [a_ncld]*len(pVal_list_N_mthyl) )
#             list_seq_pos.extend( seq_pos_list )
#             list_pVals.extend( pVal_list_N_mthyl )

#         array_out = np.transpose( np.vstack(( np.array([selected_chr_FINAL]*len(list_chr)), np.array(list_chr), np.array(list_seq_pos), np.array(list_pVals) )) )
#         print('The shape of the final array: ', array_out.shape, ' for the chr: ', selected_chr_FINAL)

#         if final_array == 0:
#             final_array = array_out
#         else:
#             final_array = np.vstack(( final_array, array_out )) 
        
#     print()
#     return final_array


# def get_values_Epinano(result_file, ref_content, col_name):

#     print()
#     print('The file is: ', result_file)

#     epinano_df = read_outputs(result_file)
#     # cheui_df = pd.read_csv(result_file[0], delimiter='\t')
#     # cheui_df = cheui_df.astype({'position':'int'})

#     selected_chr_temp = epinano_df['Ref'].unique()

#     print('For m6A: only A1 and A2 are considered. These chromosomes are overlooked: ', [x for x in selected_chr_temp if x != 'A1' and x != 'A2'])

#     selected_chr = []
#     selected_chr.extend([x for x in selected_chr_temp if x == 'A1' or x == 'A2'])

#     final_array = 0
#     for selected_chr_FINAL in natsorted(selected_chr):
#         print()
        
#         selected_chr_seq = ref_content.split('>'+selected_chr_FINAL)[1].split('>')[0]
#         print('The number of ALL nucleotides in the reference sequence: ', len(selected_chr_seq), ' for chromosome:', selected_chr_FINAL)
        
        
#         epinano_df__selected_temp = epinano_df[epinano_df['Ref'] == selected_chr_FINAL]
#         epinano_df__selected = epinano_df__selected_temp.reset_index()
#         epinano_df__selected = epinano_df__selected.drop('index', axis=1)

#         # print(epinano_df__selected)

#         list_chr = []
#         list_seq_pos = []
#         list_pVals = []
#         ncld_list = ['A', 'C', 'G', 'T']
#         for a_ncld in ncld_list:
#             # if a_chr == 'C' or a_chr == 'G' or a_chr == 'T':
#             if a_ncld == 'C' or a_ncld == 'G':
#                 continue

#             pVal_list_N_mthyl = []

#             seq_pos_list = natsorted([i for i, j in enumerate(selected_chr_seq) if j == a_ncld])
#             print('Nucleotide: ', a_ncld, ' appears ', len(seq_pos_list), ' times in the reference sequence.')

#             # As, the k-mer is a 5-mer: 

#             # print(epinano_df__selected)

#             ind_l = []
#             pos_values = []
#             for i in range(len(epinano_df__selected)):
#                 # print(epinano_df__selected['#Kmer'].values[i][2])
#                 if epinano_df__selected['#Kmer'].values[i][2] == a_ncld:
#                     x = epinano_df__selected['Window'].values[i]
#                     # print('x:', x)
#                     if x.startswith('-'):
#                         # continue
#                         pos = int((int(x.split('-')[1]) + int(x.split('-')[2]))/2)-int(x.split('-')[1])
#                     else: 
#                         pos = int((int(x.split('-')[1]) - int(x.split('-')[0]))/2)+int(x.split('-')[0])
#                     pos_values.append( pos )
#                     ind_l.append(i)
#                     # print( pos )
#                     # print(pos_values)
#             epinano_df__chr = epinano_df__selected.loc[ind_l, :]

#             # pos_values = [ int((int(epinano_df__selected['Window'].values[i].split('-')[1]) - int(epinano_df__selected['Window'].values[i].split('-')[0]))/2)+int(epinano_df__selected['Window'].values[i].split('-')[0]) for i in range(len(epinano_df__selected)) if epinano_df__selected['#Kmer'].values[i][2] == a_ncld ]
#             # print('pos_values: ', (pos_values), ' for ncld: ', a_ncld)

#             # pos_values.append(5)
#             ## Just check that these postions are correct .... 
#             [print('\n\n******** The nucleotide ', a_ncld, ' does not appear at the postion ', x, ' in the reference sequence for chr: ', selected_chr_FINAL,' as suggested by the CHEUI output. INVESTIGATE... ******\n\n') for x in pos_values if selected_chr_seq[x] != a_ncld]


#             epinano_df__chr['position'] = pos_values
#             epinano_df__chr = epinano_df__chr.drop_duplicates()
#             epinano_df__chr = epinano_df__chr.sort_values(by=['position'])

#             # print('epinano_df__chr: ')
#             # print(epinano_df__chr)

#             print('Nucleotide: ', a_ncld, ' appears ', len(epinano_df__chr), ' times in the CHEUI output 9-MERS.')

#             # print( 'pos_list: ', seq_pos_list )
#             # print( 'xPore_df__chr: ', pd.to_numeric(epinano_df__chr['position'].values) )
#             pos_DIFF = [x for x in seq_pos_list if x not in pd.to_numeric(epinano_df__chr['position'].values) ]

#             if len(pos_DIFF) > 0:
#                 print('The positions where ', a_ncld, ' appears in the reference sequence, but not in the CHEUI output: ', [x for x in pos_DIFF] )
#                 print('Total mis-match: ', len(pos_DIFF))
            
#             [pVal_list_N_mthyl.insert(x-1, pVal_prob(1)) for x in pos_DIFF]
            
#             for a_pos in seq_pos_list:
#                 for i in range(len(epinano_df__chr)):
                    
#                     if a_pos == int(epinano_df__chr['position'].values[i]) :
#                         pVal_list_N_mthyl.append(  epinano_df__chr[col_name].values[i] ) 

#             list_chr.extend( [a_ncld]*len(pVal_list_N_mthyl) )
#             list_seq_pos.extend( seq_pos_list )
#             list_pVals.extend( pVal_list_N_mthyl )

#         array_out = np.transpose( np.vstack(( np.array([selected_chr_FINAL]*len(list_chr)), np.array(list_chr), np.array(list_seq_pos), np.array(list_pVals) )) )
#         print('The shape of the final array: ', array_out.shape, ' for the chr: ', selected_chr_FINAL)

#         if final_array == 0:
#             final_array = array_out
#         else:
#             final_array = np.vstack(( final_array, array_out )) 
        
#     print()
#     return final_array



def epinano_extract(df__chr, kmer_col_name, n_mer, a_ncld, pos_col_name):
    ind_l = []
    pos_values = []
    for i in range(len(df__chr)):
        # print(epinano_df__selected['#Kmer'].values[i][2])
        if df__chr[kmer_col_name].values[i][int((n_mer-1)/2)] == a_ncld:
            x = df__chr['Window'].values[i]
            # print('x:', x)
            if x.startswith('-'):
                # continue
                pos = int((int(x.split('-')[1]) + int(x.split('-')[2]))/2)-int(x.split('-')[1])
            else: 
                pos = int((int(x.split('-')[1]) - int(x.split('-')[0]))/2)+int(x.split('-')[0])
            pos_values.append( pos )
            ind_l.append(i)
            # print( pos )
            # print(pos_values)
    df__ncld = df__chr.loc[ind_l, :]


    df__ncld[pos_col_name] = pos_values
    df__ncld = df__ncld.drop_duplicates()
    df__ncld = df__ncld.sort_values(by=[pos_col_name])

    return df__ncld, pos_values




def m6Anet_extract(df__chr, kmer_col_name, n_mer, a_ncld, pos_diff_val, pos_col_name):
    ## Get the dataframe for the selected nucleotide 
    ind_l = []
    for i in range(len(df__chr)):
        if df__chr._get_value(i, kmer_col_name)[int((n_mer-1)/2)] == a_ncld:
            ind_l.append(i)
    df__ncld = df__chr.loc[ind_l, :]

    df__ncld = df__ncld.drop_duplicates()
    pos_values = [x+pos_diff_val for x in df__ncld[pos_col_name].values]
    df__ncld[pos_col_name] = pos_values
    df__ncld = df__ncld.sort_values(by=[pos_col_name])

    return df__ncld, pos_values






def get_values_tools(tool_name, result_file, ref_content, chr_col_name, pos_col_name, kmer_col_name, prob_col_name, pos_diff_val, n_mer):

    print()
    print('The file is: ', result_file)

    if tool_name == 'CHEUI':
        df = pd.read_csv(result_file[0], delimiter='\t')
    else:
        df = read_outputs(result_file)

    if pos_col_name in df:
        df = df.astype({pos_col_name:'int'}) 
    # cheui_df = cheui_df.astype({'position':'int'})
    # print(df)

    selected_chr_temp = df[chr_col_name].unique()

    print('For m6A: only A1 and A2 are considered. These chromosomes are overlooked: ', [x for x in selected_chr_temp if x != 'A1' and x != 'A2'])

    selected_chr = []
    selected_chr.extend([x for x in selected_chr_temp if x == 'A1' or x == 'A2'])

    final_array = 0
    for selected_chr_FINAL in natsorted(selected_chr):
        print()
        
        selected_chr_seq = ref_content.split('>'+selected_chr_FINAL)[1].split('>')[0]
        print('The number of ALL nucleotides in the reference sequence: ', len(selected_chr_seq), ' for chromosome:', selected_chr_FINAL)
        
        ## Keep the dataframe with the selected chr 
        df__selected_temp = df[df[chr_col_name] == selected_chr_FINAL]
        df__chr = df__selected_temp.reset_index()
        df__chr = df__chr.drop('index', axis=1)


        list_chr = []
        list_seq_pos = []
        list_pVals = []
        ncld_list = ['A', 'C', 'G', 'T']
        for a_ncld in ncld_list:
            # if a_chr == 'C' or a_chr == 'G' or a_chr == 'T':
            if a_ncld == 'C' or a_ncld == 'G':
                continue

            pVal_list_N_mthyl = []

            seq_pos_list = natsorted([i for i, j in enumerate(selected_chr_seq) if j == a_ncld])
            print('Nucleotide: ', a_ncld, ' appears ', len(seq_pos_list), ' times in the reference sequence.')

            if tool_name == 'm6Anet' or tool_name == 'CHEUI': 
                df__ncld, pos_values = m6Anet_extract(df__chr, kmer_col_name, n_mer, a_ncld, pos_diff_val, pos_col_name)
            elif tool_name == 'Epinano':
                df__ncld, pos_values = epinano_extract(df__chr, kmer_col_name, n_mer, a_ncld, pos_col_name)

            print('Nucleotide: ', a_ncld, ' appears ', len(df__ncld), ' times in the ', tool_name, ' output ', n_mer, '-MERS.')

            output_pos_list = pd.to_numeric(df__ncld[pos_col_name].values)
            # print( 'pos_list in the seq file: ', seq_pos_list )
            # print( 'pos_list in the output file: ', output_pos_list )
            pos_DIFF = [x for x in seq_pos_list if x not in output_pos_list ]

            if len(pos_DIFF) > 0:
                print('The positions where ', a_ncld, ' appears in the reference sequence, but not in the CHEUI output: ', [x for x in pos_DIFF] )
                print('Total mis-match: ', len(pos_DIFF))
            
            [print('\n\n******** The nucleotide ', a_ncld, ' does not appear at the postion ', x, ' in the reference sequence for chr: ', selected_chr_FINAL,' as suggested by the CHEUI output. INVESTIGATE... ******\n\n') for x in pos_values if selected_chr_seq[x] != a_ncld]

            ## Add, the calculated probability values for the postions, which appears in the reference sequence but not in the tool output 
            # [pVal_list_N_mthyl.insert(x-1, pVal_prob(1)) for x in pos_DIFF]
            
            for a_pos in seq_pos_list:

                # df.query("team=='A' and position=='G'")["points"].values
                # df__ncld.query(pos_col_name==a_pos)[prob_col_name]

                
                if a_pos in output_pos_list:
                    for i in range(len(df__ncld)):
                        if a_pos == int(df__ncld[pos_col_name].values[i]) :
                            # print(df__ncld[prob_col_name].values[i])
                            pVal_list_N_mthyl.append(  df__ncld[prob_col_name].values[i] ) 
                        
                    # print(df__ncld[df__ncld[pos_col_name] == a_pos][prob_col_name])
                    # pVal_list_N_mthyl.append(  df__ncld.query(pos_col_name==a_pos)[prob_col_name].values ) 
                else:
                    pVal_list_N_mthyl.insert(a_pos-1, pVal_prob(1))

                # for i in range(len(df__ncld)):
                #     if a_pos == int(df__ncld[pos_col_name].values[i]) :
                #         # print(df__ncld[prob_col_name].values[i])
                #         pVal_list_N_mthyl.append(  df__ncld[prob_col_name].values[i] ) 



            list_chr.extend( [a_ncld]*len(pVal_list_N_mthyl) )
            list_seq_pos.extend( seq_pos_list )
            list_pVals.extend( pVal_list_N_mthyl )

        array_out = np.transpose( np.vstack(( np.array([selected_chr_FINAL]*len(list_chr)), np.array(list_chr), np.array(list_seq_pos), np.array(list_pVals) )) )
        print('The shape of the final array: ', array_out.shape, ' for the chr: ', selected_chr_FINAL)

        if final_array == 0:
            final_array = array_out
        else:
            final_array = np.vstack(( final_array, array_out )) 
        
    print()
    return final_array



def get_values_Tombo(out_file, ref_content):

    # pos_list = []
    # prob_list = []
    # for a_line in read_file(out_file).split('\n'):
    #     if len(a_line) >=1:
    #         if a_line[0].isdigit():
    #             # print(a_line.split(' ')[0], a_line.split(' ')[1])
    #             pos_list.append( a_line.split(' ')[0] )
    #             prob_list.append( a_line.split(' ')[1] )

    # outfile_lines = read_file(out_file).split('\n')
    # chr_index = []
    # for i in range(len(outfile_lines)):
    #     # if outfile_lines[i].startswith('variableStep chrom=') and selected_chr_FINAL in outfile_lines[i]:
    #     if outfile_lines[i].startswith('variableStep chrom='):
    #         chr_index.append(i)
    # print(chr_index)
    
    
    print('For m6A: only A1 and A2 are considered for Tombo.')
    selected_chr = ['A1', 'A2']

    final_array = 0
    for selected_chr_FINAL in natsorted(selected_chr):
        print()
        
        selected_chr_seq = ref_content.split('>'+selected_chr_FINAL)[1].split('>')[0]
        print('The number of ALL nucleotides in the reference sequence: ', len(selected_chr_seq), ' for chromosome:', selected_chr_FINAL)
        
        pos_list = []
        prob_list = []

        for a_line in read_file(out_file).split('variableStep chrom='+selected_chr_FINAL+' span=1')[-1].split('variableStep chrom=')[0].split('\n'):
            if len(a_line) >=1:
                # print(a_line)
                if a_line[0].isdigit():
                    # print(a_line.split(' ')[0], a_line.split(' ')[1])
                    pos_list.append( int(a_line.split(' ')[0]) )
                    prob_list.append( float(a_line.split(' ')[1]) )
        
        print('Chr: ', selected_chr_FINAL)
        print(pos_list)
        print(prob_list)


        list_chr = []
        list_seq_pos = []
        list_pVals = []
        ncld_list = ['A', 'C', 'G', 'T']
        for a_ncld in ncld_list:
            # if a_chr == 'C' or a_chr == 'G' or a_chr == 'T':
            if a_ncld == 'C' or a_ncld == 'G':
                continue

            pVal_list_N_mthyl = []

            seq_pos_list = natsorted([i for i, j in enumerate(selected_chr_seq) if j == a_ncld])
            print('Nucleotide: ', a_ncld, ' appears ', len(seq_pos_list), ' times in the reference sequence.')

            output_pos_list = []
            output_prob_list = []
            for i in range(len(pos_list)):
                if selected_chr_seq[pos_list[i]] == a_ncld:
                    # print('FINE')
                    output_pos_list.append(pos_list[i])
                    output_prob_list.append(prob_list[i])
                # else:
                #     print('XXXXXXXXX')
            
            print('Nucleotide: ', a_ncld, ' appears ', len(output_pos_list), ' times in the Tombo output ')

            # for a_pos_Tombo in pos_list:
            #     if selected_chr_seq[a_pos_Tombo] == a_ncld:
            #         # print('FINE')
            #         output_pos_list.append(a_pos_Tombo)
            #         output_prob_list.append()
            #     else:
            #         print('XXXXXXXXX')

            # # df__ncld = 
            # # col_names = ['chr', 'ref_seq_Nucl', 'position', 'DRUMMER_pVals__0', 'DRUMMER_pVals__1']
            # df__ncld = pd.DataFrame(columns=[''])

            # # df__ncld, pos_values = m6Anet_extract(df__chr, kmer_col_name, n_mer, a_ncld, pos_diff_val, pos_col_name)
            # print('Nucleotide: ', a_ncld, ' appears ', len(df__ncld), ' times in the ', tool_name, ' output ', n_mer, '-MERS.')

            # output_pos_list = pd.to_numeric(df__ncld[pos_col_name].values)
            # print( 'pos_list in the seq file: ', seq_pos_list )
            # print( 'pos_list in the output file: ', output_pos_list )
            pos_DIFF = [x for x in seq_pos_list if x not in output_pos_list ]

            if len(pos_DIFF) > 0:
                print('The positions where ', a_ncld, ' appears in the reference sequence, but not in the CHEUI output: ', [x for x in pos_DIFF] )
                print('Total mis-match: ', len(pos_DIFF))
            
            [print('\n\n******** The nucleotide ', a_ncld, ' does not appear at the postion ', x, ' in the reference sequence for chr: ', selected_chr_FINAL,' as suggested by the CHEUI output. INVESTIGATE... ******\n\n') for x in output_pos_list if selected_chr_seq[x] != a_ncld]

            ## Add, the calculated probability values for the postions, which appears in the reference sequence but not in the tool output 
            # [pVal_list_N_mthyl.insert(x-1, pVal_prob(1)) for x in pos_DIFF]
            
            for a_pos in seq_pos_list:

                if a_pos in output_pos_list:
                    for i in range(len(output_pos_list)):
                        if a_pos == int(output_pos_list[i]) :
                            # print(df__ncld[prob_col_name].values[i])
                            pVal_list_N_mthyl.append(  output_prob_list[i] ) 
                        
                    # print(df__ncld[df__ncld[pos_col_name] == a_pos][prob_col_name])
                    # pVal_list_N_mthyl.append(  df__ncld.query(pos_col_name==a_pos)[prob_col_name].values ) 
                else:
                    pVal_list_N_mthyl.insert(a_pos-1, pVal_prob(1))

                # for i in range(len(df__ncld)):
                #     if a_pos == int(df__ncld[pos_col_name].values[i]) :
                #         # print(df__ncld[prob_col_name].values[i])
                #         pVal_list_N_mthyl.append(  df__ncld[prob_col_name].values[i] ) 



            list_chr.extend( [a_ncld]*len(pVal_list_N_mthyl) )
            list_seq_pos.extend( seq_pos_list )
            list_pVals.extend( pVal_list_N_mthyl )

        array_out = np.transpose( np.vstack(( np.array([selected_chr_FINAL]*len(list_chr)), np.array(list_chr), np.array(list_seq_pos), np.array(list_pVals) )) )
        print('The shape of the final array: ', array_out.shape, ' for the chr: ', selected_chr_FINAL)

        if final_array == 0:
            final_array = array_out
        else:
            final_array = np.vstack(( final_array, array_out )) 
        
    print()
    return final_array







    





if __name__=='__main__' and '__file__' in globals():

    startTime = datetime.now()
    current_time = startTime.strftime("%Y/%m/%d at %H:%M:%S")
    print('\n\nThe script starting at: ' + str(current_time), ' \n\n' )

    ## This data token is consistent for one data type and all existing tools 
    base_dir = '/home/madhu/work/comp_tools/paper_6/download_tools'
    data_token = 'IVT_m6A'

    ref_file='/home/madhu/work/ref_transcriptome/IVT_seq/IVT_seq.fa'
    ref_content = read_file(ref_file)



    # ### DRUMMER
    # drummer_dir = base_dir+'/DRUMMER/RESULTS'
    # y_true_DRUMMER = []
    # y_score_DRUMMER = []

    # out_files__0 = [file for file in glob.glob("/home/madhu/work/comp_tools/paper_6/download_tools/DRUMMER/RESULTS/IVT_m6A_0/IVT_m6A_0__isoform/IVT_normalA_fast5_0-IVT_normalA_fast5_1/complete_analysis/*.txt")]
    # neg_array = file_content_array(natsorted(out_files__0))

    # out_files__1 = [file for file in glob.glob("/home/madhu/work/comp_tools/paper_6/download_tools/DRUMMER/RESULTS/IVT_m6A_1/IVT_m6A_1__isoform/IVT_normalA_fast5_1-IVT_m6A_fast5/complete_analysis/*.txt")]
    # pos_array = file_content_array(natsorted(out_files__1))

    # ## check that the first three columns from the two array are the same .... 


    # final_array_df_DRUMMER = np.append(neg_array, pos_array[:, 3].reshape(len(pos_array[:, 1]), 1), 1)
    # col_names = ['chr', 'ref_seq_Nucl', 'position', 'DRUMMER_pVals__0', 'DRUMMER_pVals__1']
    # df_ROC_DRUMMER = pd.DataFrame(columns=col_names)

    # df_ROC_DRUMMER = pd.DataFrame(final_array_df_DRUMMER, columns=col_names)
    # print('df_ROC_DRUMMER: ')
    # print( df_ROC_DRUMMER )












    # ### xPore 
    # xPore_dir = base_dir+'/xPore/RESULTS'
    # y_true_xPore = []
    # y_score_xPore = []


    # data_token_tmp = 'IVT_m6A_0'
    # wt_out_files = glob.glob( xPore_dir+'/'+data_token_tmp+'out_diffmod/diffmod.table' )
    # neg_array = get_values_xPore(wt_out_files, ref_content, 'pval_ko_vs_wt')

    # data_token_tmp = 'IVT_m6A_1'
    # ko_out_files = glob.glob( xPore_dir+'/'+data_token_tmp+'out_diffmod/diffmod.table' )
    # pos_array = get_values_xPore(ko_out_files, ref_content, 'pval_ko_vs_wt')


    # final_array_df_xPore = np.append(neg_array, pos_array[:, 3].reshape(len(pos_array[:, 1]), 1), 1)
    # col_names = ['chr', 'ref_seq_Nucl', 'position', 'xPore_pVals__0', 'xPore_pVals__1']
    # df_ROC_xPore = pd.DataFrame(columns=col_names)
    # df_ROC_xPore = pd.DataFrame(final_array_df_xPore, columns=col_names)
    # print('df_ROC_xPore: ')
    # print( df_ROC_xPore )



    # if len(df_ROC_DRUMMER) != len(df_ROC_xPore):
    #     print('There is something wrong in extracting values... CHECK ...')

    # df_ROC = df_ROC_DRUMMER.join(df_ROC_xPore['xPore_pVals__0'])
    # df_ROC = df_ROC.join(df_ROC_xPore['xPore_pVals__1'])

    # print('df_ROC ::::::    ')
    # print(df_ROC)


    # # Save the csv file .... 

    # save_CSV_file = '/home/madhu/work/comp_tools/paper_6/download_tools/ROC_outputs/roc_info_DRUMMER.csv'
    # df_ROC_DRUMMER.to_csv(save_CSV_file, index=False)

    # save_CSV_file = '/home/madhu/work/comp_tools/paper_6/download_tools/ROC_outputs/roc_info_xPore.csv'
    # df_ROC_xPore.to_csv(save_CSV_file, index=False)



    # save_CSV_file = '/home/madhu/work/comp_tools/paper_6/download_tools/ROC_outputs/roc_info.csv'
    # df_ROC.to_csv(save_CSV_file, index=False)





    #*******************************************************************************************


    # ### CHEUI
    # ## CHEUI states the position of 'CACTATAGC' as 30 and it starts at the position 31 (in my calculation) in the ref IVT seq: A1 
    # ## They consider the 9-mer 

    # cheui_dir = base_dir+'/CHEUI/RESULTS'

    # wt_out_files = glob.glob( cheui_dir+'/'+data_token+'/WT/WT_site_level_pred.txt' )
    # ko_out_files = glob.glob( cheui_dir+'/'+data_token+'/KO/KO_site_level_pred.txt' )

    # # cheui_df = pd.read_csv(wt_out_files[0], delimiter='\t')
    # # cheui_df = pd.read_csv(ko_out_files[0], delimiter='\t')

    # # print('cheui_df: ')
    # # print(cheui_df['position'])
    # # print(cheui_df.sort_values(by=['position']))
    # # print(cheui_df[cheui_df[['contig'=='A1']]])
    # # print(cheui_df[cheui_df['contig'] == 'A1'].sort_values(by=['position']))

    # # neg_array = get_values_CHEUI(wt_out_files, ref_content, 'probability')
    # # pos_array = get_values_CHEUI(ko_out_files, ref_content, 'probability')

    # neg_array = get_values_tools('CHEUI', wt_out_files, ref_content, 'contig', 'position', 'site', 'probability', 5, 9)
    # pos_array = get_values_tools('CHEUI', ko_out_files, ref_content, 'contig', 'position', 'site', 'probability', 5, 9)

    # final_array_df_CHEUI = np.append(neg_array, pos_array[:, 3].reshape(len(pos_array[:, 1]), 1), 1)
    # col_names = ['chr', 'ref_seq_Nucl', 'position', 'CHEUI_probs__0', 'CHEUI_probs__1']
    # df_ROC_CHEUI = pd.DataFrame(columns=col_names)
    # df_ROC_CHEUI = pd.DataFrame(final_array_df_CHEUI, columns=col_names)
    # print('df_ROC_CHEUI: ')
    # print( df_ROC_CHEUI )
    # # print(df_ROC_CHEUI[df_ROC_CHEUI['position'] == '35'])

    # save_CSV_file = '/home/madhu/work/comp_tools/paper_6/download_tools/ROC_outputs/roc_info_CHEUI.csv'
    # df_ROC_CHEUI.to_csv(save_CSV_file, index=False)






    # ### Epinano
    # ## The output file to be read 
    # epinano_dir = base_dir+'/EpiNano/RESULTS'
    # # out_files = glob.glob( epinano_dir+'/'+data_token+'/pred_DELTA*.csv' )
    # ko_out_files = glob.glob( epinano_dir+'/'+data_token+'/KO_'+data_token+'.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv' )
    # wt_out_files = glob.glob( epinano_dir+'/'+data_token+'/WT_'+data_token+'.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv' )

    # ## 'ProbM': modified probability; 'ProbU': unmodified probability

    # # neg_array = get_values_Epinano(wt_out_files, ref_content, 'ProbM')
    # # pos_array = get_values_Epinano(ko_out_files, ref_content, 'ProbM')

    # neg_array = get_values_tools('Epinano', wt_out_files, ref_content, 'Ref', 'position', '#Kmer', 'ProbM', 0, 5)
    # pos_array = get_values_tools('Epinano', ko_out_files, ref_content, 'Ref', 'position', '#Kmer', 'ProbM', 0, 5)


    # final_array_df_epinano = np.append(neg_array, pos_array[:, 3].reshape(len(pos_array[:, 1]), 1), 1)
    # col_names = ['chr', 'ref_seq_Nucl', 'position', 'Epinano_probs__0', 'Epinano_probs__1']
    # df_ROC_Epinano = pd.DataFrame(columns=col_names)
    # df_ROC_Epinano = pd.DataFrame(final_array_df_epinano, columns=col_names)
    # print('df_ROC_Epinano: ')
    # print( df_ROC_Epinano )
    # # print(df_ROC_CHEUI[df_ROC_CHEUI['position'] == '35'])

    # save_CSV_file = '/home/madhu/work/comp_tools/paper_6/download_tools/ROC_outputs/roc_info_Epinano.csv'
    # df_ROC_Epinano.to_csv(save_CSV_file, index=False)
    
    # ko_epinano_df = read_outputs(ko_out_files)
    # wt_epinano_df = read_outputs(wt_out_files)
    # # print(ko_epinano_df)
    # print(ko_epinano_df[ko_epinano_df['Window'] == '-1-3'])
    # print(ko_epinano_df[ko_epinano_df['Window'] == '0-4'])
    # print(ko_epinano_df[ko_epinano_df['Window'] == '1-5'])
    


    # ### m6Anet
    # m6Anet_dir = base_dir+'/m6anet/RESULTS'

    # wt_out_files = glob.glob( m6Anet_dir+'/'+data_token+'/results/WT/data.site_proba.csv' )
    # ko_out_files = glob.glob( m6Anet_dir+'/'+data_token+'/results/KO/data.site_proba.csv' )

    # # neg_array = get_values_Epinano(wt_out_files, ref_content, 'probability_modified')

    # neg_array = get_values_tools('m6Anet', wt_out_files, ref_content, 'transcript_id', 'transcript_position', 'kmer', 'probability_modified', 1, 5)
    # pos_array = get_values_tools('m6Anet', ko_out_files, ref_content, 'transcript_id', 'transcript_position', 'kmer', 'probability_modified', 1, 5)

    # final_array_df_m6Anet = np.append(neg_array, pos_array[:, 3].reshape(len(pos_array[:, 1]), 1), 1)
    # col_names = ['chr', 'ref_seq_Nucl', 'position', 'm6Anet_probs__0', 'm6Anet_probs__1']
    # df_ROC_m6Anet = pd.DataFrame(columns=col_names)
    # df_ROC_m6Anet = pd.DataFrame(final_array_df_m6Anet, columns=col_names)
    # print('df_ROC_m6Anet: ')
    # print( df_ROC_m6Anet )
    # # print(df_ROC_m6Anet[df_ROC_m6Anet['m6Anet_probs__0'] == '0.2356990426778793'])
    # # print(df_ROC_m6Anet[df_ROC_m6Anet['position'] == '76'])

    # save_CSV_file = '/home/madhu/work/comp_tools/paper_6/download_tools/ROC_outputs/roc_info_m6Anet.csv'
    # df_ROC_m6Anet.to_csv(save_CSV_file, index=False)


    # wt_m6Anet_df = read_outputs(wt_out_files)
    # ko_m6Anet_df = read_outputs(ko_out_files)
    # print(wt_m6Anet_df)
    # print(ko_m6Anet_df)
    # print(wt_m6Anet_df['probability_modified'].values)
    # print(ko_m6Anet_df['probability_modified'].values)



    ### Tombo

    tombo_dir = base_dir+'/Tombo_Output/RESULTS'

    wt_out_files = glob.glob( tombo_dir+'/'+data_token+'/IVT_m6Afinal_text__.fraction_modified_reads.m6A.plus.wig' )
    ko_out_files = glob.glob( tombo_dir+'/'+data_token+'/IVT_m6Afinal_text__.fraction_modified_reads.m6A.minus.wig' )
    
    neg_array = get_values_Tombo(wt_out_files[0], ref_content)
    pos_array = get_values_Tombo(ko_out_files[0], ref_content)

    final_array_df_Tombo = np.append(neg_array, pos_array[:, 3].reshape(len(pos_array[:, 1]), 1), 1)
    col_names = ['chr', 'ref_seq_Nucl', 'position', 'Tombo_probs__0', 'Tombo_probs__1']
    df_ROC_Tombo = pd.DataFrame(columns=col_names)
    df_ROC_Tombo = pd.DataFrame(final_array_df_Tombo, columns=col_names)
    print('df_ROC_Tombo: ')
    print( df_ROC_Tombo )
    # print(df_ROC_m6Anet[df_ROC_m6Anet['m6Anet_probs__0'] == '0.2356990426778793'])
    print(df_ROC_Tombo[df_ROC_Tombo['position'] == '186'])

    save_CSV_file = '/home/madhu/work/comp_tools/paper_6/download_tools/ROC_outputs/roc_info_Tombo.csv'
    df_ROC_Tombo.to_csv(save_CSV_file, index=False)









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

