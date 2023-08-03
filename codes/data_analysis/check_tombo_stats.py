
import os,sys
import math


sys.path.insert(0, '/home/qian/anaconda3/envs/ont3/lib/python3.6/site-packages')
#/home/qian/anaconda3/envs/ont3/lib/python3.6/site-packages/tombo/tombo_stats.py

from Bio import SeqIO
from tombo import tombo_stats
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import metrics

opn = 15

import copy

def calculate_perf( cm ):
    num_true = cm[0][0] + cm[0][1];
    num_fals = cm[1][0] + cm[1][1];

    total = num_true + num_fals
    pred_t = cm[0][0] + cm[1][0]
    pred_f = cm[0][1] + cm[1][1]

    recall = cm[0][0]/(1 if num_true==0 else num_true)
    precision = cm[0][0]/(1 if pred_t==0 else pred_t)
    acc = (cm[0][0] + cm[1][1])/(1 if total==0 else total)
    f1 = 2*recall*precision/(recall + precision if recall + precision>0 else 1)

    print("R={:3f} P=={:3f} Acc={:3f} F1={:3f}".format( recall, precision, acc, f1 ))



def get_min_max(mlist):
    cplist = sorted(copy.deepcopy(mlist))
    
    min_ = cplist[ int(len(cplist)*0.05+0.5) ]
    max_ = cplist[ int(len(cplist)*0.95+0.5) ]

    return min_, max_

def get_percentile(tlist, a_pos):
    binnum = 20;

    rrach_list = []
    for _r1 in ['A', 'G']:
        for _r2 in ['A', 'G']:
            for _h1 in ['A', 'C', 'T']:
                rrach_list.append( ''.join( [_r1, _r2, 'AC', _h1] ) );
    all_a = []
    for _r1 in ['A', 'G', 'C', 'T']:
        for _r2 in ['A', 'G', 'C', 'T']:
            for _r4 in ['A', 'G', 'C', 'T']:
                for _r5 in ['A', 'G', 'C', 'T']:
                    all_a.append( ''.join( [_r1, _r2, 'A', _r4, _r5] ) )

    
    for_auc = [[], []]

    epi_predict_files = ['Curlcake_m6a/Curlcake_m6a.plus_strand.merge.per.site.csv.modification.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv', 'Curlcake_non/Curlcake_non.plus_strand.merge.per.site.csv.modification.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv']
    ep_label_svm = ['EpinanoSVM.A', 'EpinanoSVM.RRACH']
    ep_m_num_svm = []
    pred_result = [[], []]
    for _ef in epi_predict_files:
        orignalv = []
        mod_unmod = [[0, 0], [0, 0]]
        pred_res = pd.read_csv( _ef, header=0);
        print("All", pred_res.shape)
        pred_res_r = pred_res[ pred_res['Ref']=='cc6m_2244_T7_ecorv' ]
        print('Ref', pred_res_r.shape);
        pred_res_rA    = pred_res_r[ pred_res_r['#Kmer'].isin(all_a) ] 
        print("Aal", pred_res_rA.shape)
        pred_res_rrach = pred_res_r[ pred_res_r['#Kmer'].isin(rrach_list) ]
        print('RRA', pred_res_rrach.shape)
        mod_unmod[0][0] = (pred_res_rA['prediction']=='mod').sum()
        mod_unmod[0][1] =  (pred_res_rA['prediction']=='unm').sum()
        mod_unmod[1][0] = (pred_res_rrach['prediction']=='mod').sum()
        mod_unmod[1][1] =  (pred_res_rrach['prediction']=='unm').sum()

        pred_result[0].append( mod_unmod[0] )
        pred_result[1].append( mod_unmod[1] )

        if 'm6a' in _ef:
            for_auc[0].append( ( [1 for _ in range(pred_res_rA['ProbM'].shape[0])], pred_res_rA['ProbM'].tolist() ) )
            for_auc[1].append( "EpiNano" )
        else:
            for_auc[0][-1][0].extend( [0 for _ in range(pred_res_rA['ProbM'].shape[0]) ]  )
            for_auc[0][-1][1].extend( pred_res_rA['ProbM'].tolist() )
            print("*********************len: {} {}".format(  len(for_auc[0][-1][0]), len( for_auc[0][-1][1] ) ))

        if 'm6a' in _ef:
            orignalv = pred_res_rA['ProbM'].tolist()
            print(len(orignalv), orignalv[:3])
            min_, max_ = get_min_max( orignalv )
            m1_num = [0 for _ in range(binnum)]
            for pp_ind in range(len( orignalv )):
               bini = 1-(orignalv[pp_ind]-min_)/(max_-min_)
               ind_p = int( 100*bini/(100/binnum)  +0.5)
               if ind_p<0:
                   print("test: ", min_, max_, orignalv[pp_ind], bini, ind_p, binnum)
                   ind_p = 0
               for _dp in range(ind_p, binnum):
                   m1_num[_dp] += 1;
            for _i in range(binnum):
                m1_num[_i] = m1_num[_i]/(len(a_pos));
            ep_m_num_svm.append( m1_num )
            orignalv = pred_res_rrach['ProbM'].tolist()
            print(len(orignalv), orignalv[:3])
            min_, max_ = get_min_max( orignalv )
            m1_num = [0 for _ in range(binnum)]
            for pp_ind in range(len( orignalv )):
               bini = 1-(orignalv[pp_ind]-min_)/(max_-min_)
               ind_p = int( 100*bini/(100/binnum)  +0.5)
               if ind_p<0:
                   print("test: ", min_, max_, orignalv[pp_ind], bini, ind_p, binnum)
                   ind_p = 0
               for _dp in range(ind_p, binnum):
                   m1_num[_dp] += 1;
            for _i in range(binnum):
                m1_num[_i] = m1_num[_i]/(len(a_pos));
            ep_m_num_svm.append( m1_num )
    calculate_perf(pred_result[0])
    calculate_perf(pred_result[1])

    epinano_op_files = ['Test.delta-sum_err.prediction.csv', 'Test.linear-regression.prediction.csv']
    ep_percentile = []
    ep_m_num = []
    ep_label = ['Epinano.Delta', 'Epinano.Linear']
    for _ef in epinano_op_files:
        orignalv = []
        mod_unmod = [0, 0]
        with open(_ef, 'r') as mr:
            line = mr.readline();
            while True:
                line = mr.readline();
                if not line: break;
                line = line.strip();
                lsp = line.split(',')
                if len(lsp)<1: continue;
                c_p_info = lsp[0].split();
                if (c_p_info[0] not in ['cc6m_2244_T7_ecorv']) or (c_p_info[-1] not in ['+']) or (c_p_info[-2] not in ['A', 'a']): continue
                #print( len(lsp), lsp )
                if len(lsp)<6: continue
                if lsp[-1] in ['unm']: mod_unmod[1] += 1
                else: 
                    mod_unmod[0] += 1
                    print(line)
                orignalv.append( float( "{:.2f}".format( float( lsp[-2] ) ) ) )
        print(len(orignalv), orignalv[:10], mod_unmod) 
        min_, max_ = get_min_max( orignalv )
        m1_num = [0 for _ in range(binnum)]
        for pp_ind in range(len( orignalv )):
            bini = 1-(orignalv[pp_ind]-min_)/(max_-min_)
            ind_p = int( 100*bini/(100/binnum)  +0.5)
            if ind_p<0:
                print("test: ", min_, max_, orignalv[pp_ind], bini, ind_p, binnum)
                ind_p = 0
            for _dp in range(ind_p, binnum):
                m1_num[_dp] += 1;
        for _i in range(binnum):
            m1_num[_i] = m1_num[_i]/(len(a_pos));
        ep_m_num.append( m1_num )

    min_ = None;
    max_ = None;
    dlist = []
    for pp_ind in range(len(tlist)):
        if tlist[pp_ind][1] not in a_pos: continue;
        dlist.append( float("{:.2f}".format(tlist[pp_ind][0])) )
        if min_==None or min_>tlist[pp_ind][0]:
            min_ = tlist[pp_ind][0]
        if max_==None or max_<tlist[pp_ind][0]:
            max_ = tlist[pp_ind][0]
    min_, max_ = get_min_max( dlist )
    dlist = sorted(dlist)
    print('\n', dlist[:opn], dlist[-opn:])

    binnum = 20;
    m_num = [0 for _ in range(binnum)]
    for pp_ind in range(len(tlist)):
        if tlist[pp_ind][1] not in a_pos: continue;
        bini = 1-(tlist[pp_ind][0]-min_)/(max_-min_)
        ind_p = int( 100*bini/(100/binnum)  +0.5)
        if ind_p<0:
           print("test.t: ", min_, max_, tlist[pp_ind], bini, ind_p, binnum)
           ind_p = 0
        for _dp in range(ind_p, binnum):
            m_num[_dp] += 1;
    for _i in range(binnum):
        m_num[_i] = m_num[_i]/(len(a_pos));
   
    m_percentile = [_i/100.0 for _i in range(0, 100, 100//binnum)]

    print( m_percentile )
    print( m_num )

    plt.rcParams["figure.figsize"] = [8, 6]
    plt.rcParams["figure.autolayout"] = True

    plt.plot( m_percentile, ep_m_num[0], 'b-.', label=ep_label[0]);
    plt.plot( m_percentile, ep_m_num[1], 'g-.', label=ep_label[1]);
    #plt.plot( m_percentile, ep_m_num_svm[0], 'r-.', label=ep_label_svm[0]);
    #plt.plot( m_percentile, ep_m_num_svm[1], 'k-.', label=ep_label_svm[1]);


    ax1 = plt.subplot()
    plt.plot( m_percentile, m_num, 'k-', label="Tombo")
    #plt.plot( m_percentile_new, m_num_new, label="NEW")

    otherpred = ['log.TS_finetune_pred.l3h1024lr10000ep1.20222.log', 'log.TS_finetune_pred.l3h1024lr10000p1.ep9.26820.log', 'log.TS_finetune_pred.l3h1024lr10000ep2.675674.log']
    #otherpred = ['log.TS_finetune_pred.l3h1024lr10000ep2.675674.log', 'log.TS_finetune_pred.l3h1024lr10000p1.ep9.26820.log', 'log.TS_finetune_pred.l3h1024lr10000ep1.20222.log']
    colormst = ['r', 'b', 'g']
    lables = ['fine-tune', '1% fine-tune', 'fine-tune 2']
    min_ = 0
    max_ = 1
    for ofid_idx in range(len(otherpred)-1):
        ofid = otherpred[ ofid_idx ]
        m_num = [0 for _ in range(binnum)]

        for_t_roc = [[], []]

        with open(ofid, 'r') as mr:
            start_str = 'Pred: cc6m_2244_T7_ecorv:+:'
            perclist = []
            perclist2 = []
            read_perf = [[0,0],[0,0]]
            site_perf = [[0,0],[0,0]]
            site_perf_h = [[0,0],[0,0]]
            while True:
                line = mr.readline();
                if not line: break;
                line = line.strip();
                if len(line)<len(start_str): continue;
                if not line[:len(start_str)]==start_str: continue;

                lsp = line[len(start_str):].split()
                if int(lsp[0]) not in a_pos: print("Error A not in: {}".format(lsp[0]))
                perclist.append( float(lsp[1])/float(int(lsp[1]) + int(lsp[2])) )
                perclist2.append( float(lsp[3])/float(int(lsp[3]) + int(lsp[4])) )

                read_perf[0][0] += float(lsp[1])
                read_perf[0][1] += float(lsp[2])
                read_perf[1][0] += float(lsp[3])
                read_perf[1][1] += float(lsp[4])
                if float(lsp[1])/float(int(lsp[1]) + int(lsp[2]))<0.25:
                    site_perf[0][1] += 1
                else: site_perf[0][0] += 1
                if float(lsp[3])/float(int(lsp[3]) + int(lsp[4]))<0.25:
                    site_perf[1][1] += 1
                else: site_perf[1][0] += 1

                for_t_roc[1].append( float(lsp[1])/float(int(lsp[1]) + int(lsp[2])) )
                for_t_roc[0].append( 1 )
                for_t_roc[1].append( float(lsp[3])/float(int(lsp[3]) + int(lsp[4])) )
                for_t_roc[0].append( 0 )

                if float(lsp[1])/float(int(lsp[1]) + int(lsp[2]))<0.5:
                    site_perf_h[0][1] += 1
                else: site_perf_h[0][0] += 1
                if float(lsp[3])/float(int(lsp[3]) + int(lsp[4]))<0.5:
                    site_perf_h[1][1] += 1
                else: site_perf_h[1][0] += 1

                ind_p = int( (1-perclist[-1])*100/(100/binnum)  +0.5)
                for _dp in range(ind_p, binnum):
                  m_num[_dp] += 1;
            print(len(for_auc), len(for_auc[0]), len( for_auc[1] ) )
            for_auc[0].append( (for_t_roc[0], for_t_roc[1]) )
            for_auc[1].append( lables[ofid_idx] )
            print("*********************len: {} {}".format(  len(for_auc[0][-1][0]), len( for_auc[0][-1][1] ) ))
            print( ofid )
            print("site_perf: {}".format(site_perf))
            calculate_perf(site_perf)
            print("site_perfH: {}".format(site_perf_h))
            calculate_perf(site_perf_h)
            print("read_perf: {}".format(read_perf))
            calculate_perf( read_perf )
            read_perf[0][0] /= 1e6 
            read_perf[0][1] /= 1e6
            read_perf[1][0] /= 1e6
            read_perf[1][1] /= 1e6
            for _i in range(2):
                for _j in range(2):
                    read_perf[_i][_j] = float( "{:.2f}".format(read_perf[_i][_j]) )
            print("read_perf: {}".format(read_perf))

            for _i in range(binnum):
                m_num[_i] = m_num[_i]/(len(a_pos));
            #plt.plot( m_percentile, m_num, colormst[ofid_idx]+'-', label=ofid[len('log.TS_finetune_pred.'):-len('.log')])
            min_ = None;
            max_ = None;
            rtlist = []
            for pp_ind in range(len(perclist)):
                ratio = math.log( perclist[pp_ind]/perclist2[pp_ind] )
                rtlist.append( float("{:.2f}".format(ratio) ) )
                if min_==None or min_>ratio:
                    min_ = ratio
                if max_==None or max_<ratio:
                    max_ = ratio
            #print(rtlist[:5])
            min_, max_ = get_min_max( rtlist )
            #print(rtlist[:5])
            rtlist = sorted(rtlist)
            m_num = [0 for _ in range(binnum)]
            for pp_ind in range(len(perclist)):
                ratio = math.log( perclist[pp_ind]/perclist2[pp_ind] )
                bini = 1-(ratio-min_)/(max_-min_)
                ind_p = int( 100*bini/(100/binnum)  +0.5)
                if ind_p<0:
                   print("testR"+lables[ofid_idx], min_, max_, perclist[pp_ind], bini, ind_p, binnum)
                   ind_p = 0
                for _dp in range(ind_p, binnum):
                    m_num[_dp] += 1;
            for _i in range(binnum):
                m_num[_i] = m_num[_i]/(len(a_pos));
            #plt.plot( m_percentile, m_num, colormst[ofid_idx]+'-.', label=ofid[len('log.TS_finetune_pred.'):-len('.log')]+'.ratio')
            min_, max_ = get_min_max( perclist )
            m_num = [0 for _ in range(binnum)]
            for pp_ind in range(len(perclist)):
                bini = 1-(perclist[pp_ind]-min_)/(max_-min_)
                ind_p = int( 100*bini/(100/binnum)  +0.5)
                if ind_p<0:
                   print("testP"+lables[ofid_idx], min_, max_, perclist[pp_ind], bini, ind_p, binnum)
                   ind_p = 0
                for _dp in range(ind_p, binnum):
                    m_num[_dp] += 1;
            for _i in range(binnum):
                m_num[_i] = m_num[_i]/(len(a_pos));
            plt.plot( m_percentile, m_num, colormst[ofid_idx]+'--', label=lables[ofid_idx]); #ofid[len('log.TS_finetune_pred.'):-len('.log')]+'.perc')

            perclist = sorted(perclist)
            print('\n', [float("{:.2f}".format(tpl)) for tpl in perclist[:opn]], [float("{:.2f}".format(tpl)) for tpl in perclist[-opn:]] )
            print('\n', rtlist[:opn], rtlist[-opn:])

    mlegend = plt.legend(fontsize=10, ncol=2)
    ax1.set_xlim([0, 1]);
    plt.xlabel('Score rank (higher to lower from left to right)')
    plt.ylabel('Percentage')
    #plt.show()
    plt.savefig('check_tombo_stats.png', dpi=600)

    plt.figure()
    plt.rcParams["figure.figsize"] = [8, 6]
    plt.rcParams["figure.autolayout"] = True
    for _aui in range(len(for_auc[0])):
        fpr, tpr, _ = metrics.roc_curve( for_auc[0][_aui][0] , for_auc[0][_aui][1]  )
        auc = metrics.roc_auc_score(for_auc[0][_aui][0] , for_auc[0][_aui][1] )
        plt.plot(fpr,tpr, colormst[_aui]+'--', label=("{}, auc={:.3f}".format( for_auc[1][_aui] , auc)) )

    plt.xlabel('1-Specificity (FPR)')
    plt.ylabel('Sensitivity(Recall)')
    mlegend = plt.legend(fontsize=10, ncol=2)
    plt.savefig('check_tombo_stats_auc.png', dpi=600)



if __name__ == '__main__':
    m_g = SeqIO.index(sys.argv[1], 'fasta')

    #epinano_op = ['Test.delta-sum_err.prediction.csv', 'Test.linear-regression.prediction.csv']

    if sys.argv[3]=='l':
        per_read_stats = tombo_stats.LevelStats(sys.argv[2])
    else:
        per_read_stats = tombo_stats.ModelStats(sys.argv[2])


    for seq_ in m_g:
        if seq_ not in ['cc6m_2244_T7_ecorv']: continue
        a_pos = []
        other_pos = []
        #print(m_g[seq_].seq)   
        for _pind in range( len(m_g[seq_].seq) ):
            #print ( m_g[seq_].seq[ _pind ], m_g[seq_].seq[ _pind ] in ['A', 'a'] ) 
            if m_g[seq_].seq[ _pind ] in ['A', 'a']:
                a_pos.append( _pind )
        print( a_pos )
        per_pos_stats = per_read_stats.get_reg_stats(seq_, '+',1,len(m_g[seq_].seq))
        #print( per_pos_stats )
        for pp_ind in range(len(per_pos_stats)):
            if pp_ind<10:
                print( per_pos_stats[pp_ind] )
            else: break;
        if sys.argv[3]=='l':
            get_percentile(per_pos_stats , a_pos)



