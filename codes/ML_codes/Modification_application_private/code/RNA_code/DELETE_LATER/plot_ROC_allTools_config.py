



CHR_SELECTED = ['C1', 'C2']
NCLD_SELECTED = ['C', 'G']





# import os,sys
# # import time, datetime
# sys.path.insert(0, '/home/madhu/work/codes')
# from main import check_endWith, check_create_dir, check_env, read_file, write_file, count_lines, save_file, load_file

# # check_env('torch')

# # import TS_finetune_pred_MP
# from datetime import datetime
# import glob

# import warnings
# # warnings.filterwarnings("ignore", category=FutureWarning)
# warnings.filterwarnings("ignore")

# from sklearn.metrics import roc_auc_score, roc_curve
# import matplotlib.pyplot as plt
# import numpy as np
# import pandas as pd
# from tqdm import tqdm
# import pickle

# from natsort import natsorted

# '''
# This script draws the ROC curves for the other existing tools. 
# It reads the output from those tools. Tools are: 

# 1. Epinano
# 2. m6Anet
# 3. Tombo 
# 4. DRUMMER
# 5. xPore
# 6. CHEUI
# 7. Our Tool
# '''

# def read_outputs(csv_files):
#     if len(csv_files) != 1:
#         print('******** There are multiple csv files. INVESTIGATE ... *********')
#     else:
#         csv_file = csv_files[0]
    
#     # print('The final output file is: ', csv_file)
#     df = pd.read_csv(csv_file)
#     return df




# # def draw_ROC(y_true, y_score, save_ROC_results, saveFig_title):
# # values is a list containing: [y_true, y_score, tool]
# def calc_AUC__draw_ROC(values, saveFig_title): 

#     fig, ax = plt.subplots(figsize=(20, 20))

#     list_c = ['r', 'b', 'm', 'y', 'g', 'c', 'k']
#     list_m = ['*', 'x', 'o', '.', '*', 'x', '*']
#     list_s = ['-', '--', '-.', ':', ' ', '', '-']
#     # list_s = ['-', '--', '-.', ':', ' ', '', 'solid', 'dashed', 'dashdot', 'dotted']

#     if len(values) != len(list_c) or len(values) != len(list_c):
#         print('***** SOMETHING is WRONG, INVESTIGATE ... ******')
#         sys.exit(0)
    
#     # for a_value in values:
#     for i in range(len(values)): 
#         [y_true, y_score, tool] = values[i]
#         fpr, tpr, thresholds = roc_curve(y_true, y_score)
#         AUC_score = roc_auc_score(y_true, y_score)
#         print('\n\nThe AUC is :', round(AUC_score, 4), ' for the TOOL: ', tool)
#         # plt.plot(fpr, tpr, list_c[i], marker=list_m[i], linestyle=list_s[i], linewidth=2, label='{}: (AUC = {:.4f})'.format(tool, AUC_score))
#         plt.plot(fpr, tpr, list_c[i], linewidth=3, label='{}: (AUC = {:.4f})'.format(tool, AUC_score))

#     plt.legend(loc='lower right', fontsize=28)
#     plt.plot([0, 1], [0, 1], color='navy', linestyle='--')
#     plt.xlim([-0.02, 1.02])
#     plt.ylim([-0.02, 1.02])
#     plt.xlabel('Specificity', fontsize=25)
#     plt.ylabel('Sensitivity', fontsize=25)
#         # plt.xlabel('False Positive Rate', fontsize=20)
#         # plt.ylabel('True Positive Rate', fontsize=20)
    
#     n_split = 11
#     xTick_list = []
#     for n in np.linspace(0, 1, n_split):
#         xTick_list.append(str(int(n*100))+'%')
#     # reversing the list
#     new_xTick_list = []
#     for i in xTick_list:
#         new_xTick_list.insert(0, i)
#     plt.xticks(np.linspace(0, 1, n_split), new_xTick_list, fontsize=22)
#     yTick_list = []
#     for n in np.linspace(0, 1, n_split):
#         yTick_list.append(str(int(n*100))+'%')
#     plt.yticks(np.linspace(0, 1, n_split), yTick_list, fontsize=22)
#     plt.grid(color='y', linewidth=0.5)
    
#     # plt.title('Mean ROC curve for TB Index Score', fontsize=35)
#     plt.show()
#     plt.savefig(saveFig_title)
#     print('The ROC figure has been saved at: ', saveFig_title)
#     plt.close('all')
    
#     # save_file(save_ROC_results, (fpr, tpr, thresholds, AUC_score))

#     # return AUC_score



# def file_content_array(tool_name, out_files, ref_content, chr_col_name, pos_col_name, kmer_col_name, prob_col_name, pos_diff_val, n_mer, if_pVal):
#     print('out_files: ', out_files)
#     for i in range(len(out_files)):
#         array_out = get_values_tools(tool_name, out_files[i], ref_content, chr_col_name, pos_col_name, kmer_col_name, prob_col_name, pos_diff_val, n_mer, if_pVal)
#         if array_out != 0:
#             if i == 0:
#                 array = array_out
#             else:
#                 array = np.vstack(( array, array_out )) 
#     return array





# def pVal_prob(p_val):
#     # -np.log( np.float_(drummer_df__chr[col_name].values[i]) if np.float_(drummer_df__chr[col_name].values[i])>10**(-300) else 10**(-300))
#     prob = -np.log( np.float_(p_val) if np.float_(p_val)>10**(-300) else 10**(-300))
#     return prob



# def extract_epinano(df__chr, kmer_col_name, n_mer, a_ncld, pos_col_name):
#     ind_l = []
#     pos_values = []
#     for i in range(len(df__chr)):
#         # print(epinano_df__selected['#Kmer'].values[i][2])
#         if df__chr[kmer_col_name].values[i][int((n_mer-1)/2)] == a_ncld:
#             x = df__chr['Window'].values[i]
#             # print('x:', x)
#             if x.startswith('-'):
#                 # continue
#                 pos = int((int(x.split('-')[1]) + int(x.split('-')[2]))/2)-int(x.split('-')[1])
#             else: 
#                 pos = int((int(x.split('-')[1]) - int(x.split('-')[0]))/2)+int(x.split('-')[0])
#             pos_values.append( pos )
#             ind_l.append(i)
#             # print( pos )
#             # print(pos_values)
#     df__ncld = df__chr.loc[ind_l, :]


#     df__ncld[pos_col_name] = pos_values
#     df__ncld = df__ncld.drop_duplicates()
#     df__ncld = df__ncld.sort_values(by=[pos_col_name])

#     return df__ncld, pos_values




# def extract_tools(df__chr, kmer_col_name, n_mer, a_ncld, pos_diff_val, pos_col_name, if_DRUMMER):
#     ## Get the dataframe for the selected nucleotide 
#     ind_l = []

#     for i in range(len(df__chr)):
#         if if_DRUMMER == 0:
#             if df__chr._get_value(i, kmer_col_name)[int((n_mer-1)/2)] == a_ncld:
#                 ind_l.append(i)
#         elif if_DRUMMER == 1:
#         # for i in range(len(df__chr)):
#             if df__chr._get_value(i, kmer_col_name) == a_ncld:
#                 ind_l.append(i)
#     df__ncld = df__chr.loc[ind_l, :]

#     df__ncld = df__ncld.drop_duplicates()
#     pos_values = [x+pos_diff_val for x in df__ncld[pos_col_name].values]
#     df__ncld[pos_col_name] = pos_values
#     df__ncld = df__ncld.sort_values(by=[pos_col_name])

#     return df__ncld, pos_values



# def check_pos_diff(seq_pos_list, output_pos_list, a_ncld, tool_name, selected_chr_FINAL, pos_values, selected_chr_seq):
#     # print( 'pos_list in the seq file: ', seq_pos_list )
#     # print( 'pos_list in the output file: ', output_pos_list )
#     pos_DIFF = [x for x in seq_pos_list if x not in output_pos_list ]
#     if len(pos_DIFF) > 0:
#         # print('The positions where ', a_ncld, ' appears in the reference sequence, but not in the '+tool_name+' output: ', [x for x in pos_DIFF] )
#         print('Total mis-match: ', len(pos_DIFF))
#         if len(output_pos_list)+len(pos_DIFF) == len(seq_pos_list):
#             print('The number of appearnaces in ref seq - tool output is OK. ')
#         else:
#             print('\n\n******** PLEASE CHECK THIS ******\n\n')

#     [print('\n\n******** The nucleotide ', a_ncld, ' does not appear at the postion ', x, ' in the reference sequence for chr: ', selected_chr_FINAL,' as suggested by the CHEUI output. INVESTIGATE... ******\n\n') for x in pos_values if selected_chr_seq[x] != a_ncld]



# ## if_pVal = 0, if probability is given
# ## if_pVal = 1, if p-value is given
# def get_values_tools(tool_name, result_file, ref_content, chr_col_name, pos_col_name, kmer_col_name, prob_col_name, pos_diff_val, n_mer, if_pVal):

#     print()
#     print('The file is: ', result_file)

#     if tool_name == 'CHEUI':
#         df = pd.read_csv(result_file[0], delimiter='\t')
#     elif tool_name == 'DRUMMER':
#         file_content = []
#         for a_line in read_file(result_file).split('\n'):
#             if not len(a_line.split('\t')) == 1:
#                 file_content.append(a_line.split('\t'))
#         df = pd.DataFrame(data=np.array(file_content))
#         df.columns = df.iloc[0]
#         df = df.reset_index(drop=True)
#         df = df[1:]
#     else:
#         df = read_outputs(result_file)

#     if pos_col_name in df:
#         df = df.astype({pos_col_name:'int'}) 

#     # print(df)

#     if tool_name == 'DRUMMER':
#          selected_chr = [result_file.split('/')[-1].split('.')[0]]
#     else:
#         selected_chr = CHR_SELECTED
    
#     print('selected_chr: ', selected_chr)

#     final_array = 0
#     for selected_chr_FINAL in natsorted(selected_chr):
#         print()
        
#         selected_chr_seq = ref_content.split('>'+selected_chr_FINAL)[1].split('>')[0]
#         print('The number of ALL nucleotides in the reference sequence: ', len(selected_chr_seq), ' for chromosome:', selected_chr_FINAL)
        
#         ## Keep the dataframe with the selected chr 
#         df__selected_temp = df[df[chr_col_name] == selected_chr_FINAL]
#         df__chr = df__selected_temp.reset_index()
#         df__chr = df__chr.drop('index', axis=1)

#         # print(df__chr)
#         # print(df__chr[df__chr['pos_ctrl'] == 22])

#         list_chr = []
#         list_seq_pos = []
#         list_pVals = []
#         ncld_list = ['A', 'C', 'G', 'T']
#         for a_ncld in ncld_list:
#             # if a_chr == 'C' or a_chr == 'G' or a_chr == 'T':
#             # if a_ncld == 'C' or a_ncld == 'G':
#             if not a_ncld in NCLD_SELECTED:
#                 continue

#             pVal_list = []

#             seq_pos_list = natsorted([i for i, j in enumerate(selected_chr_seq) if j == a_ncld])
#             print('Nucleotide: ', a_ncld, ' appears ', len(seq_pos_list), ' times in the reference sequence.')

#             if tool_name == 'm6Anet' or tool_name == 'CHEUI' or tool_name == 'xPore': 
#                 df__ncld, pos_values = extract_tools(df__chr, kmer_col_name, n_mer, a_ncld, pos_diff_val, pos_col_name, 0)
#             elif tool_name == 'DRUMMER':
#                 df__ncld, pos_values = extract_tools(df__chr, kmer_col_name, n_mer, a_ncld, pos_diff_val, pos_col_name, 1)
#             elif tool_name == 'Epinano':
#                 df__ncld, pos_values = extract_epinano(df__chr, kmer_col_name, n_mer, a_ncld, pos_col_name)
            

#             print('Nucleotide: ', a_ncld, ' appears ', len(df__ncld), ' times in the ', tool_name, ' output ', n_mer, '-MERS.')

#             output_pos_list = pd.to_numeric(df__ncld[pos_col_name].values)

#             check_pos_diff(seq_pos_list, output_pos_list, a_ncld, tool_name, selected_chr_FINAL, pos_values, selected_chr_seq)

#             ## Add, the calculated probability values for the postions, which appears in the reference sequence but not in the tool output 
#             # [pVal_list.insert(x-1, pVal_prob(1)) for x in pos_DIFF]
            
#             for a_pos in seq_pos_list:
                
#                 if a_pos in output_pos_list:
#                     for i in range(len(df__ncld)):
#                         if a_pos == int(df__ncld[pos_col_name].values[i]) :
#                             # print(df__ncld[prob_col_name].values[i])
#                             if if_pVal == 0:
#                                 pVal_list.append(  df__ncld[prob_col_name].values[i] ) 
#                             elif if_pVal == 1:
#                                 pVal_list.append(  pVal_prob(df__ncld[prob_col_name].values[i]) ) 
#                 else:
#                     pVal_list.insert(a_pos-1, pVal_prob(1))


#             list_chr.extend( [a_ncld]*len(pVal_list) )
#             list_seq_pos.extend( seq_pos_list )
#             list_pVals.extend( pVal_list )

#         array_out = np.transpose( np.vstack(( np.array([selected_chr_FINAL]*len(list_chr)), np.array(list_chr), np.array(list_seq_pos), np.array(list_pVals) )) )
#         print('The shape of the final array: ', array_out.shape, ' for the chr: ', selected_chr_FINAL)

#         if final_array == 0:
#             final_array = array_out
#         else:
#             final_array = np.vstack(( final_array, array_out )) 
        
#     print()
#     return final_array



# def get_values_Tombo(tool_name, out_file, ref_content):
    
#     # print('For m6A: only A1 and A2 are considered for Tombo.')
#     # selected_chr = ['A1', 'A2']
#     # selected_chr = CHR_SELECTED

#     final_array = 0
#     for selected_chr_FINAL in natsorted(CHR_SELECTED):
#         print()
        
#         selected_chr_seq = ref_content.split('>'+selected_chr_FINAL)[1].split('>')[0]
#         print('The number of ALL nucleotides in the reference sequence: ', len(selected_chr_seq), ' for chromosome:', selected_chr_FINAL)
        
#         pos_list = []
#         prob_list = []

#         for a_line in read_file(out_file).split('variableStep chrom='+selected_chr_FINAL+' span=1')[-1].split('variableStep chrom=')[0].split('\n'):
#             if len(a_line) >=1:
#                 # print(a_line)
#                 if a_line[0].isdigit():
#                     # print(a_line.split(' ')[0], a_line.split(' ')[1])
#                     pos_list.append( int(a_line.split(' ')[0]) )
#                     prob_list.append( float(a_line.split(' ')[1]) )
        
#         print('Chr: ', selected_chr_FINAL)
#         print(pos_list)
#         print(prob_list)


#         list_chr = []
#         list_seq_pos = []
#         list_pVals = []
#         ncld_list = ['A', 'C', 'G', 'T']
#         for a_ncld in ncld_list:
#             # if a_chr == 'C' or a_chr == 'G' or a_chr == 'T':
#             # if a_ncld == 'C' or a_ncld == 'G':
#             if not a_ncld in NCLD_SELECTED:
#                 continue

#             pVal_list = []

#             seq_pos_list = natsorted([i for i, j in enumerate(selected_chr_seq) if j == a_ncld])
#             print('Nucleotide: ', a_ncld, ' appears ', len(seq_pos_list), ' times in the reference sequence.')

#             output_pos_list = []
#             output_prob_list = []
#             for i in range(len(pos_list)):
#                 if selected_chr_seq[pos_list[i]] == a_ncld:
#                     output_pos_list.append(pos_list[i])
#                     output_prob_list.append(prob_list[i])

            
#             print('Nucleotide: ', a_ncld, ' appears ', len(output_pos_list), ' times in the ', tool_name, ' output ')

#             check_pos_diff(seq_pos_list, output_pos_list, a_ncld, tool_name, selected_chr_FINAL, output_pos_list, selected_chr_seq)

#             ## Add, the calculated probability values for the postions, which appears in the reference sequence but not in the tool output 
#             # [pVal_list.insert(x-1, pVal_prob(1)) for x in pos_DIFF]
            
#             for a_pos in seq_pos_list:

#                 if a_pos in output_pos_list:
#                     for i in range(len(output_pos_list)):
#                         if a_pos == int(output_pos_list[i]) :
#                             # print(df__ncld[prob_col_name].values[i])
#                             pVal_list.append(  output_prob_list[i] ) 
#                 else:
#                     pVal_list.insert(a_pos-1, pVal_prob(1))

#             list_chr.extend( [a_ncld]*len(pVal_list) )
#             list_seq_pos.extend( seq_pos_list )
#             list_pVals.extend( pVal_list )

#         array_out = np.transpose( np.vstack(( np.array([selected_chr_FINAL]*len(list_chr)), np.array(list_chr), np.array(list_seq_pos), np.array(list_pVals) )) )
#         print('The shape of the final array: ', array_out.shape, ' for the chr: ', selected_chr_FINAL)

#         if final_array == 0:
#             final_array = array_out
#         else:
#             final_array = np.vstack(( final_array, array_out )) 
        
#     print()
#     return final_array



# def get_values_ourtool(tool_name, result_file, ref_content):

#     print()

#     # print('For m6A: only A1 and A2 are considered for Tombo.')
#     # selected_chr = ['A1', 'A2']
#     # selected_chr = CHR_SELECTED

#     final_array = 0
#     # for selected_chr_FINAL in natsorted(selected_chr):
#     for selected_chr_FINAL in natsorted(CHR_SELECTED):
#         print()
        
#         selected_chr_seq = ref_content.split('>'+selected_chr_FINAL)[1].split('>')[0]
#         print('The number of ALL nucleotides in the reference sequence: ', len(selected_chr_seq), ' for chromosome:', selected_chr_FINAL)

#         pos_list = []
#         strand_list = []
#         prob_list_0 = []
#         prob_list_1 = []
#         for a_part in (result_file.split('\n')):
#             if a_part.startswith('Pred: '+selected_chr_FINAL):
#                 info_part = a_part.split(' ')
#                 if len(info_part) == 6:
#                     # print('a_part: ', a_part)
#                     strand = info_part[1].split(':')[1]
#                     pos = int(info_part[1].split(':')[2])
#                     conf_mat = info_part[2:len(info_part)]
#                     # pos_list.append(pos+1)
#                     # strand_list.append(strand)
                    
#                     if (int(conf_mat[0])+int(conf_mat[1])) != 0 and (int(conf_mat[2])+int(conf_mat[3])) != 0:
#                         pos_list.append(pos+1)
#                         strand_list.append(strand)
#                         # y_true_ours.append(0)
#                         prob_list_0.append(int(conf_mat[1])/(int(conf_mat[0])+int(conf_mat[1])))
#                         # y_true_ours.append(1)
#                         prob_list_1.append(int(conf_mat[3])/(int(conf_mat[2])+int(conf_mat[3])))
#                     # else:
#                     #     print('a_part: ', a_part)
        
#         # print('pos_list: ', len(pos_list))
#         # print('prob_list_0: ', len(prob_list_0))
#         # print('prob_list_1: ', len(prob_list_1))
#         df_ourTool = pd.DataFrame(np.transpose( np.vstack(( np.array([selected_chr_FINAL]*len(pos_list)), np.array(pos_list), np.array(strand_list), np.array(prob_list_0), np.array(prob_list_1) )) ), columns=['Chr', 'pos', 'strand', 'prob_list_0', 'prob_list_1'])
#         df_ourTool['pos'] = pd.to_numeric(df_ourTool['pos'].values)
#         df_ourTool_temp = df_ourTool.sort_values('pos')
#         df_ourTool = df_ourTool_temp.reset_index()
#         df_ourTool = df_ourTool.drop('index', axis=1)

#         list_chr = []
#         list_seq_pos = []
#         list_pVals_0 = []
#         list_pVals_1 = []
#         ncld_list = ['A', 'C', 'G', 'T']
#         for a_ncld in ncld_list:
#             # if a_chr == 'C' or a_chr == 'G' or a_chr == 'T':
#             # if a_ncld == 'C' or a_ncld == 'G':
#             if not a_ncld in NCLD_SELECTED:
#                 continue
            
#             pVal_list_mthyl = []
#             pVal_list_N_mthyl = []

#             seq_pos_list = natsorted([i for i, j in enumerate(selected_chr_seq) if j == a_ncld])
#             print('Nucleotide: ', a_ncld, ' appears ', len(seq_pos_list), ' times in the reference sequence.')

#             output_pos_list = []
#             output_prob_list_0 = []
#             output_prob_list_1 = []
#             for i in range(len(df_ourTool['pos'])):
#                 if (df_ourTool._get_value(i, 'strand') == '+' and a_ncld == NCLD_SELECTED[0]) or (df_ourTool._get_value(i, 'strand') == '-' and a_ncld == NCLD_SELECTED[1]):
#                 # if (df_ourTool._get_value(i, 'strand') == '+' and a_ncld == 'A') or (df_ourTool._get_value(i, 'strand') == '-' and a_ncld == 'T'):
#                     output_pos_list.append(df_ourTool._get_value(i, 'pos'))
#                     output_prob_list_0.append(df_ourTool._get_value(i, 'prob_list_0'))
#                     output_prob_list_1.append(df_ourTool._get_value(i, 'prob_list_1'))
            
#             print('Nucleotide: ', a_ncld, ' appears ', len(output_pos_list), ' times in the ', tool_name, ' output ')

#             check_pos_diff(seq_pos_list, output_pos_list, a_ncld, tool_name, selected_chr_FINAL, output_pos_list, selected_chr_seq)

#             ## Add, the calculated probability values for the postions, which appears in the reference sequence but not in the tool output 
#             # [pVal_list_N_mthyl.insert(x-1, pVal_prob(1)) for x in pos_DIFF]
            
#             for a_pos in seq_pos_list:

#                 if a_pos in output_pos_list:
#                     for i in range(len(output_pos_list)):
#                         if a_pos == int(output_pos_list[i]) :
#                             # print(df__ncld[prob_col_name].values[i])
#                             pVal_list_N_mthyl.append(  output_prob_list_0[i] ) 
#                             pVal_list_mthyl.append(  output_prob_list_1[i] ) 
#                 else:
#                     pVal_list_N_mthyl.insert(a_pos-1, pVal_prob(1))
#                     pVal_list_mthyl.insert(a_pos-1, pVal_prob(1))

#             list_chr.extend( [a_ncld]*len(pVal_list_N_mthyl) )
#             list_seq_pos.extend( seq_pos_list )
#             list_pVals_0.extend( pVal_list_N_mthyl )
#             list_pVals_1.extend( pVal_list_mthyl )

#         array_out = np.transpose( np.vstack(( np.array([selected_chr_FINAL]*len(list_chr)), np.array(list_chr), np.array(list_seq_pos), np.array(list_pVals_0), np.array(list_pVals_1) )) )
#         print('The shape of the final array: ', array_out.shape, ' for the chr: ', selected_chr_FINAL)

#         if final_array == 0:
#             final_array = array_out
#         else:
#             final_array = np.vstack(( final_array, array_out )) 

#     print()
#     return final_array


# def get_true_scores(df_ROC, tool_name):
#     y_true = []
#     y_score = []

#     col_0 = tool_name+'_probs__0'
#     col_1 = tool_name+'_probs__1'

#     y_true.extend([0]*len(df_ROC))
#     y_score.extend(df_ROC[col_0])

#     y_true.extend([1]*len(df_ROC))
#     y_score.extend(df_ROC[col_1])

#     return y_true, y_score



# def get_final_array(tool_name, wt_file, ko_file, ref_content, chr_col_name, pos_col_name, kmer_col_name, prob_col_name, pos_diff_val, n_mer, if_pVal):
#     if tool_name == 'Tombo':
#         neg_array = get_values_Tombo(tool_name, wt_file[0], ref_content)
#         pos_array = get_values_Tombo(tool_name, ko_file[0], ref_content)
#     elif tool_name == 'OUR':
#         log_file_content = read_file(wt_file)
#         # ref_gene_content = read_file(ko_file)
#         final_array_df = get_values_ourtool(tool_name, log_file_content, ref_content)
#     elif tool_name == 'DRUMMER':
#         neg_array = file_content_array(tool_name, natsorted(wt_file), ref_content, chr_col_name, pos_col_name, kmer_col_name, prob_col_name, pos_diff_val, n_mer, if_pVal)
#         pos_array = file_content_array(tool_name, natsorted(ko_file), ref_content, chr_col_name, pos_col_name, kmer_col_name, prob_col_name, pos_diff_val, n_mer, if_pVal)
#     else:
#         neg_array = get_values_tools(tool_name, wt_file, ref_content, chr_col_name, pos_col_name, kmer_col_name, prob_col_name, pos_diff_val, n_mer, if_pVal)
#         pos_array = get_values_tools(tool_name, ko_file, ref_content, chr_col_name, pos_col_name, kmer_col_name, prob_col_name, pos_diff_val, n_mer, if_pVal)
    

#     if tool_name != 'OUR':
#         final_array_df = np.append(neg_array, pos_array[:, 3].reshape(len(pos_array[:, 1]), 1), 1)
    
#     col_names = ['chr', 'ref_seq_Nucl', 'position', tool_name+'_probs__0', tool_name+'_probs__1']
#     df_ROC_tool = pd.DataFrame(final_array_df, columns=col_names)
#     print('The final df for the tool: ', tool_name)
#     print( df_ROC_tool )
    
#     ## Save the csv file .... 
#     save_CSV_file = '/home/madhu/work/comp_tools/paper_6/download_tools/ROC_outputs/roc_info_'+tool_name+'.csv'
#     df_ROC_tool.to_csv(save_CSV_file, index=False)

#     return df_ROC_tool



# def get_DRUMMER_files(file_path):
#     DRUMMER_wt = []
#     for file in glob.glob(file_path):
#         for a_chr in CHR_SELECTED:
#             if a_chr+'.' in file:
#                 DRUMMER_wt.append(file)
#     return DRUMMER_wt


# def get_AUC(y_true, y_score):
#     AUC_score = roc_auc_score(y_true, y_score)
#     print('\n\n**** AUC: ', AUC_score, '****\n\n')
#     return AUC_score



# if __name__=='__main__' and '__file__' in globals():

#     startTime = datetime.now()
#     current_time = startTime.strftime("%Y/%m/%d at %H:%M:%S")
#     print('\n\nThe script starting at: ' + str(current_time), ' \n\n' )

#     ## This data token is consistent for one data type and all existing tools 
#     base_dir = '/home/madhu/work/comp_tools/paper_6/download_tools'



#     #####################################################################################
#     '''
#     This one is about IVT m6A data
#     '''
#     data_token = 'IVT_m6A'
#     ref_file='/home/madhu/work/ref_transcriptome/IVT_seq/IVT_seq.fa'
#     ref_content = read_file(ref_file)
#     CHR_SELECTED = ['A1', 'A2']
#     NCLD_SELECTED = ['A', 'T']

#     drummer_dir = base_dir+'/DRUMMER/RESULTS'
#     DRUMMER_wt = get_DRUMMER_files( drummer_dir+"/IVT_m6A_0/IVT_m6A_0__isoform/IVT_normalA_fast5_0-IVT_normalA_fast5_1/complete_analysis/*.txt" )
#     DRUMMER_ko = get_DRUMMER_files( drummer_dir+"/IVT_m6A_1/IVT_m6A_1__isoform/IVT_normalA_fast5_1-IVT_m6A_fast5/complete_analysis/*.txt" )

#     xPore_dir = base_dir+'/xPore/RESULTS'
#     data_token_tmp = 'IVT_m6A_0'
#     xPore_wt = glob.glob( xPore_dir+'/'+data_token_tmp+'out_diffmod/diffmod.table' )
#     data_token_tmp = 'IVT_m6A_1'
#     xPore_ko = glob.glob( xPore_dir+'/'+data_token_tmp+'out_diffmod/diffmod.table' )

#     our_log_file = '/home/madhu/Logs/predict__train_p3_BOTH_pred_p5_m6A_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.240025.pt__A.log'    
#     # our_log_file = '/home/madhu/Logs/predict__train_p3_BOTH_pred_p5_m6A_A_BilstmMean.layers3.hs1024.F.lr1500.b512.p1000.GPU.ep1.240025.pt__A.log'

#     #####################################################################################

#     #####################################################################################

#     '''
#     This one is about m6A data in paper 3 
#     '''
#     # data_token = 'm6A_p3__49_50'
#     # ref_file='/home/madhu/work/ref_transcriptome/OTHERS/synthetic_construct/Best_Practices_dRNAseq_analysis/reference_fasta/cc.fasta'
#     # ref_content = read_file(ref_file)

#     # NCLD_SELECTED = ['A', 'T']
#     # CHR_SELECTED = []
#     # selected_chr = []
#     # for a_part in (ref_content.split('\n')):
#     #     if a_part.startswith('>'):
#     #         # selected_chr.append( a_part.split('>')[1] )
#     #         CHR_SELECTED.append( a_part.split('>')[1] )
    
#     # print('CHR_SELECTED: ', CHR_SELECTED)

#     # drummer_dir = base_dir+'/DRUMMER/RESULTS'
#     # DRUMMER_wt = [file for file in glob.glob(drummer_dir+"/m6A_p3__49_51/m6A_p3__49_51__isoform/GSM3528749-GSM3528751/complete_analysis/*.txt")]
#     # DRUMMER_ko = [file for file in glob.glob(drummer_dir+"/m6A_p3__49_50/m6A_p3__49_50__isoform/GSM3528749-GSM3528750/complete_analysis/*.txt")]
#     # # DRUMMER_ko = [file for file in glob.glob(drummer_dir+"/m6A_p3__51_52/m6A_p3__51_52__isoform/GSM3528751-GSM3528752/complete_analysis/*.txt")]

#     # xPore_dir = base_dir+'/xPore/RESULTS'
#     # data_token_tmp = 'm6A_p3__49_51'
#     # xPore_wt = glob.glob( xPore_dir+'/'+data_token_tmp+'out_diffmod/diffmod.table' )
#     # data_token_tmp = 'm6A_p3__49_50'
#     # xPore_ko = glob.glob( xPore_dir+'/'+data_token_tmp+'out_diffmod/diffmod.table' )

#     # # our_log_file = '/home/madhu/Logs/predict__pred_p3_train_p5_m6A_ALL_A_BilstmMean.layers3.hs512.F.lr1000.b256.p1000.GPU.ep1.270255.pt__A.log'
#     # our_log_file = '/home/madhu/Logs/predict__pred_p3_train_p5_m6A_ALL_A_BilstmMean.layers3.hs1024.F.lr500.b64.p1000.GPU.ep1.310658.pt__A.log'

#     #####################################################################################

    
#     cheui_dir = base_dir+'/CHEUI/RESULTS'
#     cheui_wt = glob.glob( cheui_dir+'/'+data_token+'/WT/WT_site_level_pred.txt' )
#     cheui_ko = glob.glob( cheui_dir+'/'+data_token+'/KO/KO_site_level_pred.txt' )

#     epinano_dir = base_dir+'/EpiNano/RESULTS'
#     epinano_wt = glob.glob( epinano_dir+'/'+data_token+'/WT_'+data_token+'.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv' )
#     epinano_ko = glob.glob( epinano_dir+'/'+data_token+'/KO_'+data_token+'.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv' )

#     m6Anet_dir = base_dir+'/m6anet/RESULTS'
#     m6Anet_wt = glob.glob( m6Anet_dir+'/'+data_token+'/results/WT/data.site_proba.csv' )
#     m6Anet_ko = glob.glob( m6Anet_dir+'/'+data_token+'/results/KO/data.site_proba.csv' )

#     tombo_dir = base_dir+'/Tombo_Output/RESULTS'
#     tombo_wt = glob.glob( tombo_dir+'/'+data_token+'/'+data_token+'final_text__.fraction_modified_reads.m6A.plus.wig' )
#     tombo_ko = glob.glob( tombo_dir+'/'+data_token+'/'+data_token+'final_text__.fraction_modified_reads.m6A.minus.wig' )





#     ### DRUMMER
#     df_ROC_DRUMMER = get_final_array('DRUMMER', DRUMMER_wt, DRUMMER_ko, ref_content, 'chr_ctrl', 'pos_ctrl', 'ref_ctrl', 'p_values_OR_adj', 0, 0, 1)

#     ### xPore 
#     df_ROC_xPore = get_final_array('xPore', xPore_wt, xPore_ko, ref_content, 'id', 'position', 'kmer', 'pval_ko_vs_wt', 1, 5, 1)

#     ### CHEUI
#     ## CHEUI states the position of 'CACTATAGC' as 30 and it starts at the position 31 (in my calculation) in the ref IVT seq: A1 
#     ## They consider the 9-mer 
#     df_ROC_CHEUI = get_final_array('CHEUI', cheui_wt, cheui_ko, ref_content, 'contig', 'position', 'site', 'probability', 5, 9, 0)

#     ### Epinano
#     df_ROC_Epinano = get_final_array('Epinano', epinano_wt, epinano_ko, ref_content, 'Ref', 'position', '#Kmer', 'ProbM', 0, 5, 0)

#     ### m6Anet
#     df_ROC_m6Anet = get_final_array('m6Anet', m6Anet_wt, m6Anet_ko, ref_content, 'transcript_id', 'transcript_position', 'kmer', 'probability_modified', 1, 5, 0)

#     ### Tombo
#     df_ROC_Tombo = get_final_array('Tombo', tombo_wt, tombo_ko, ref_content, 0,0,0,0,0,0,0)

#     ### OUR DNN TOOL
#     df_ROC_OUR = get_final_array('OUR', our_log_file, 0, ref_content, 0,0,0,0,0,0,0)






#     df_ROC = df_ROC_DRUMMER.join(df_ROC_xPore['xPore_probs__0']).join(df_ROC_xPore['xPore_probs__1']).join(df_ROC_CHEUI['CHEUI_probs__0']).join(df_ROC_CHEUI['CHEUI_probs__1']).join(df_ROC_Epinano['Epinano_probs__0']).join(df_ROC_Epinano['Epinano_probs__1']).join(df_ROC_m6Anet['m6Anet_probs__0']).join(df_ROC_m6Anet['m6Anet_probs__1']).join(df_ROC_Tombo['Tombo_probs__0']).join(df_ROC_Tombo['Tombo_probs__1']).join(df_ROC_OUR['OUR_probs__0']).join(df_ROC_OUR['OUR_probs__1'])
#     # df_ROC = df_ROC_DRUMMER.join(df_ROC_xPore['xPore_probs__0']).join(df_ROC_xPore['xPore_probs__1']).join(df_ROC_CHEUI['CHEUI_probs__0']).join(df_ROC_CHEUI['CHEUI_probs__1']).join(df_ROC_Epinano['Epinano_probs__0']).join(df_ROC_Epinano['Epinano_probs__1']).join(df_ROC_m6Anet['m6Anet_probs__0']).join(df_ROC_m6Anet['m6Anet_probs__1']).join(df_ROC_Tombo['Tombo_probs__0']).join(df_ROC_Tombo['Tombo_probs__1'])
#     # df_ROC = df_ROC_DRUMMER

#     print('df_ROC ::::::    ')
#     print(df_ROC)

#     save_CSV_file = '/home/madhu/work/comp_tools/paper_6/download_tools/ROC_outputs/roc_info_FINAL_'+data_token+'.csv'
#     df_ROC.to_csv(save_CSV_file, index=False)

    




#     # ## Load the CSV file: 
#     df_ROC = pd.read_csv(save_CSV_file)

#     # # print(df_ROC)



#     y_true_DRUMMER, y_score_DRUMMER = get_true_scores(df_ROC, 'DRUMMER')
#     y_true_xPore, y_score_xPore = get_true_scores(df_ROC, 'xPore')
#     y_true_CHEUI, y_score_CHEUI = get_true_scores(df_ROC, 'CHEUI')
#     y_true_epinano, y_score_epinano = get_true_scores(df_ROC, 'Epinano')
#     y_true_m6Anet, y_score_m6Anet = get_true_scores(df_ROC, 'm6Anet')
#     y_true_Tombo, y_score_Tombo = get_true_scores(df_ROC, 'Tombo')
#     y_true_ours, y_score_ours = get_true_scores(df_ROC, 'OUR')





#     # get_AUC(y_true_DRUMMER, y_score_DRUMMER)
#     # get_AUC(y_true_xPore, y_score_xPore)
#     # get_AUC(y_true_CHEUI, y_score_CHEUI)
#     # get_AUC(y_true_epinano, y_score_epinano)
#     # get_AUC(y_true_m6Anet, y_score_m6Anet)



#     saveFig_title = base_dir + '/ROC_'+data_token+'.png'
#     # values is a list containing: [y_true, y_score, tool]
#     values = [y_true_ours, y_score_ours, 'Our_DNN_tool'], [y_true_epinano, y_score_epinano, 'Epinano'], [y_true_m6Anet, y_score_m6Anet, 'm6Anet'], [y_true_Tombo, y_score_Tombo, 'Tombo'], [y_true_DRUMMER, y_score_DRUMMER, 'DRUMMER'], [y_true_xPore, y_score_xPore, 'xPore'], [y_true_CHEUI, y_score_CHEUI, 'CHEUI']
#     calc_AUC__draw_ROC(values, saveFig_title)



#     executionTime = (datetime.now() - startTime)
#     current_time = datetime.now().strftime("%Y/%m/%d at %H:%M:%S")
#     print('\n\nThe script completed at: ' + str(current_time))
#     print('Execution time: ' + str(executionTime), ' \n\n')


