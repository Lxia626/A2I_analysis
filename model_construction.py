#! /usr/bin/env python

import numpy as np
import pandas as pd
import sys, os
import warnings
warnings.filterwarnings('ignore')
import math
from time import time

# ML models
import joblib
# from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
# from sklearn.grid_search import GridsearchCV
# from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
# from sklearn.linear_model import LinearRegression, LogisticRegression
# import sklearn.naive_bayes as sk_bayes
# from sklearn.tree import DecisionTreeClassifier
# from sklearn.svm import SVC
# from sklearn.neural_network import MLPClassifier
# from sklearn.neighbors import KNeighborsClassifier
# ML parameters
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.metrics import mean_squared_error, accuracy_score, precision_score, recall_score, f1_score, roc_auc_score
from sklearn import metrics
from sklearn.metrics import roc_curve, auc
# plots
import matplotlib.pyplot as plt
# import seaborn as sns

def arg(argv):

    import argparse
    parser = argparse.ArgumentParser(description="Calculate the site on each ALU repeat by global alignment according to RepBase data", 
                                     formatter_class=argparse.RawTextHelpFormatter, 
                                     epilog='Note:\n\tOptional models are HistGradientBoostingClassifier, RandomForestClassifier, LinearRegression, LogisticRegression, DecisionTreeClassifier, KNeighborsClassifier, multinomial_BAYES, bernoulli_BAYES, gaussian_BAYES, SVM')
    parser.add_argument('-i', '--input', help="Input data including features of all candidate sites", type=str, required=True)
    parser.add_argument('-P', '--prefix', help="Prefix for output files", type=str, required=True)
    parser.add_argument('--model_selection', help="Whether to select machine learning model. [FUNCTION, False]", action='store_true', default=False)
    parser.add_argument('--model', help="Model to use [FUNCTION, NA; if to use trained model, then supply a file with suffix as '.pkl']", type=str, default='NA')
    parser.add_argument('--bcf_out', help="Output of BCFtools in table format with mismatch number calculated [Optinal supplied for '--model' only]", type=str, default=None)
    parser.add_argument('--rand_time', help="Time to randomly select training sets [for '--model_selection' Function only, 50]", type=int, default=50)
    parser.add_argument('-p', help="Maximum number of subprocesses [for '--model_selection' Function only, 1]", type=int, default=1)
    parser.add_argument('-f', '--filter_features', help="Whether to filter features based on feature_importances_ [for '--model' Function only when training models, False]", action='store_true', default=False)
    parser.add_argument('-F', '--selected_features', help="Selected features for downstream prediction [for '--model' Function only when loading models to predict, False]", type=str, default=None)
    parser.add_argument('--set_divide', help="How to divide the dataset (training:validation:test) [for '--model' Function only, '7,2,1']", default='7,2,1', type=str)
    # parser.add_argument('-R', help="Reference sequences [for '--model' Function only]", type=str)
    # parser.add_argument('--gene_annotation', help="Annotation file for genes [for '--model' Function only]", type=str)
    parser.add_argument('--seq_error_p', help="Sequencing error [for '--model' Function only]", type=float, default=0.05)
    parser.add_argument('-ks', '--known_sites', help="File storing known editing sites. If not supply, then enter with 'none'", type=str, default='/public/home/Songlab/songyl/database/editing_sites/Merged/human/human_all_merged_editing-sites_add-strand')
    parser.add_argument('-snp', '--snp_sites', help="File storing SNP sites. If not supply, then enter with 'none'", type=str, default='/public/home/Songlab/xial/Work/database/hg19/SNP_all_dbSNP_1000Genomes_UWash.txt')
    parser.add_argument('--test_size', help="Ratio of test set", type=float, default=0.3)

    args = parser.parse_args() 
    sys.stderr.write('args: '+str(args)+'\n')

    if args.known_sites == 'none':
        ks_list = None; sys.stderr.write('Attention: No known sites were supplied.\n')
    else:
        ks_list = []
        for ks in open(args.known_sites, 'r'):
            ks_list.append(ks.split('\t')[0]+'_'+ks.split('\t')[1])
    if args.snp_sites == 'none':
        snp_list = None; sys.stderr.write('Attention: No SNP sites were supplied.\n')
    else:
        snp_list = []
        for snp in open(args.snp_sites):
            snp_list.append(snp.split('\t')[0]+'_'+snp.strip().split('\t')[1])
    # if args.seq_error_p < 1: seq_error_p = args.seq_error_p
    # else: seq_error_p = 10**((-1)*(round(args.seq_error_p/10, 1)))

    if args.model_selection:

        # data preparation
        # inputFile = '/public/home/Songlab/xial/Work/Nanopore_RNA_editing/data/own_data/N2324962_80-1411948157_2024-02-01/results/HEK293T/own_data_HEK293T_nano_150_7/04_callSite/own_data_HEK293T_nano_a2T3_2.all.ALU.candsite.txt'
        data = pd.read_csv(args.input, sep='\t', header=0, na_values=['nan'])
        data['CHROM'] = data['CHROM'].astype('str'); data['POS'] = data['POS'].astype('str')
        data['ID'] = data['CHROM'] + '_' + data['POS']
        data['KNOWN'] = 2
        data.loc[data['DEL_NUM_1'] == 0, 'DEL_LEN_1'] = 0
        data.loc[data['DEL_SITE_NUM_1'] == 0, 'DEL_SITE_RATIO_1'] = 0
        data.loc[data['SKIP_1'] == 1, 'DEL_SITE_RATIO_1'] = 0
        data.loc[data['INS_NUM_1'] == 0, 'INS_LEN_1'] = 0
        data.loc[data['DEL_NUM_mid'] == 0, 'DEL_LEN_mid'] = 0
        data.loc[data['DEL_SITE_NUM_mid'] == 0, 'DEL_SITE_RATIO_mid'] = 0
        data.loc[data['SKIP_mid'] == 1, 'DEL_SITE_RATIO_mid'] = 0
        data.loc[data['INS_NUM_mid'] == 0, 'INS_LEN_mid'] = 0
        data.loc[data['DEL_NUM_3'] == 0, 'DEL_LEN_3'] = 0
        data.loc[data['DEL_SITE_NUM_3'] == 0, 'DEL_SITE_RATIO_3'] = 0
        data.loc[data['SKIP_3'] == 1, 'DEL_SITE_RATIO_3'] = 0
        data.loc[data['INS_NUM_3'] == 0, 'INS_LEN_3'] = 0
        data['SPLICE_1'] = data['SPLICE_1'] + 20
        data['SPLICE_mid'] = data['SPLICE_mid'] + 20
        data['SPLICE_3'] = data['SPLICE_3'] + 20
        data.loc[data['ID'].isin(ks_list), 'KNOWN'] = 1
        data.loc[data['ID'].isin(snp_list), 'KNOWN'] = 0
        sys.stderr.write('#All data: %s; true_num: %s; false_num: %s; ' % (data.shape, data[data['KNOWN'] == 1].shape[0], data[data['KNOWN'] == 0].shape[0]))
        data = data[data['KNOWN'] != 2]
        data = data.drop(columns=['SEQ_ERROR_P_1', 'SEQ_ERROR_P_mid', 'SEQ_ERROR_P_3'])
        sys.stderr.write('#Filtered data: %s\n' % (data.shape[0]))
        sys.stderr.write('Features for analysing: '+str(data.iloc[:,5:57].columns.values)+'\n')
        # total_num = data.shape[0]; true_num = data[data['KNOWN'] == 1].shape[0]; false_num = data[data['KNOWN'] == 0].shape[0]

        if args.p == 1:
            selectModel(data.copy(), args.prefix+'.txt', test_size=args.test_size, rand_state=[1, args.rand_time+1])
        if args.p > 1:
            from multiprocessing import Pool
            pool = Pool(processes=args.p)
            splice = math.ceil(args.rand_time / args.p)
            sys.stderr.write('rand_time: %s; processes: %s\n' % (args.rand_time, args.p))
            rand_state = list(range(args.rand_time))
            for i in range(args.p-1):
                if i*splice+1 > len(rand_state): continue
                elif (i+1)*splice+1 > len(rand_state): _rand_stat = [rand_state[i*splice+1], args.rand_time+1]
                else: _rand_stat = [rand_state[i*splice+1], rand_state[(i+1)*splice+1]]
                sys.stderr.write('_rand_stat for process %s: %s\n' % (i, _rand_stat))
                pool.apply_async(selectModel, (data.copy(), args.prefix+'.'+str(i)+'.txt', args.test_size, _rand_stat, ))
            if (args.p-1)*splice < args.rand_time:
                _rand_stat_final = [rand_state[(args.p-1)*splice+1], args.rand_time+1]
                pool.apply_async(selectModel, (data.copy(), args.prefix+'.'+str(args.p)+'.txt', args.test_size, _rand_stat_final, ))
                sys.stderr.write('_rand_stat_final: '+str(_rand_stat_final)+'\n')
            # sys.exit()
            pool.close()
            pool.join()
            os.system('cat %s | grep -v \'Model\' | sed \'1i Model\tRound\tTrain_score\tTest_score\tAc_score\tPr_score\tRc_score\tF1_score\tROC_AUC_score\tRunningTime\' > %s' % (args.prefix+'.*.txt', args.out))
            os.system('rm '+args.prefix+'.*.txt')
    
    elif args.model:

        if args.model == 'NA': 
            sys.stderr.write('Model to use should be specified. Options are listed in help document.\n'); sys.exit()
        if args.model.split('.')[-1] != 'pkl' and args.selected_features: 
            sys.stderr.write('"-F, --selected_features" should not be applied when loading models for prediction.\n'); sys.exit()
        # from pysam import FastaFile
        # if not args.R: sys.stderr.write('Reference sequences should be supplied in "--model" mode.\n'); sys.exit()
        # ref_in = FastaFile(args.R)
        # from python_analysis import gtf_chrom_trans, base_pair
        # if not args.gene_annotation: 
        #     if not args.R: sys.stderr.write('Annotation file should be supplied in "--model" mode.\n'); sys.exit()
        # chrom_info_ori_anno = gtf_chrom_trans(args.gene_annotation)
        # chrom_anno = {}
        # for _c in chrom_info_ori_anno:
        #     if not _c.startswith('chr'): chrom = 'chr'+_c
        #     else: chrom = _c
        #     chrom_anno[chrom] = sorted(chrom_info_ori_anno[_c], key=lambda x:x[1][0]) # sort according to the end of each annotation
        # cand_site_info = {}
        # for _x in open(args.cand_sites, 'r'):
        #     if _x.startswith('CHROM'): continue
        #     cand_site_info.setdefault(_x.split('\t')[0], []).append((int(_x.split('\t')[1]), _x.split('\t')[2], _x.split('\t')[3], _x.split('\t')[4], float(_x.split('\t')[5]), int(_x.split('\t')[35]))) # (pos, reads_num, ref, alt, edit_level, depth)
        # if args.seq_error_p < 1: seq_error_p = args.seq_error_p
        # else: seq_error_p = 10**((-1)*(round(args.seq_error_p/10, 1)))
        if args.bcf_out:
            bcfout_dict = {}
            for line in open(args.bcf_out):
                if line.startswith('CHROM'): continue
                info = line.strip().split('\t')
                bcfout_dict[info[0]+"_"+info[1]] = info[4:] # {'SITE_ID': [DP, VDB, RPBZ, MQBZ, NMBZ, SCBZ, MQ0F, DP4, MQ, MN]}
        else: bcfout_dict = None

        modelConstructAndPredict(args.input, args.prefix, bcfout=bcfout_dict, ks_list=ks_list, snp_list=snp_list, test_size=args.test_size, model2Use=args.model, filter_features=args.filter_features, set_divide=args.set_divide, selected_features_str=args.selected_features)

def selectModel(data_all, outFile, test_size=0.3, rand_state=None):

    from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
    from sklearn.linear_model import LinearRegression, LogisticRegression
    import sklearn.naive_bayes as sk_bayes
    from sklearn.tree import DecisionTreeClassifier
    from sklearn.svm import SVC
    from sklearn.neural_network import MLPClassifier
    from sklearn.neighbors import KNeighborsClassifier

    out_file = open(outFile, 'w')
    out_file.write('Model\tRound\tTrain_score\tTest_score\tAc_score\tPr_score\tRc_score\tF1_score\tROC_AUC_score\tRunningTime\n')

    # data preparation
    data_target = data_all['KNOWN']
    data_features = data_all.iloc[:,5:57]
    # sys.exit()
    for _i in range(rand_state[0], rand_state[-1]):
        sys.stderr.write('Current random time: %s\n' % (_i))
        X_train, X_test, y_train, y_test = train_test_split(data_features, data_target, test_size=test_size, random_state=_i)
        START_TIME = time()
        sys.stderr.write('Model: RandomForestClassifier\n')
        train_score, test_score, ac_score, pr_score, rc_score, f_score, roc_auc_score = RandomForestClassifier_func(X_train, X_test, y_train, y_test)
        END_1 = time()
        out_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('RandomForestClassifier', _i, train_score, test_score, ac_score, pr_score, rc_score, f_score, roc_auc_score, END_1-START_TIME))
        sys.stderr.write('Model: HistGradientBoostingClassifier\n')
        train_score, test_score, ac_score, pr_score, rc_score, f_score, roc_auc_score = HistGradientBoostingClassifier_func(X_train, X_test, y_train, y_test)
        END_2 = time()
        out_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('HistGradientBoostingClassifier', _i, train_score, test_score, ac_score, pr_score, rc_score, f_score, roc_auc_score, END_2-END_1))
        sys.stderr.write('Model: RandomForestClassifier\n')
        train_score, test_score, ac_score, pr_score, rc_score, f_score, roc_auc_score = RandomForestClassifier_func(X_train, X_test, y_train, y_test)
        END_2 = time()
        out_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('RandomForestClassifier', _i, train_score, test_score, ac_score, pr_score, rc_score, f_score, roc_auc_score, END_2-END_1))
        sys.stderr.write('Model: LinearRegression\n')
        train_score, test_score, ac_score, pr_score, rc_score, f_score, roc_auc_score = LinearRegression_func(X_train, X_test, y_train, y_test)
        END_3 = time()
        out_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('LinearRegression', _i, train_score, test_score, ac_score, pr_score, rc_score, f_score, roc_auc_score, END_3-END_2))
        sys.stderr.write('Model: LogisticRegression\n')
        train_score, test_score, ac_score, pr_score, rc_score, f_score, roc_auc_score = LogisticRegression_func(X_train, X_test, y_train, y_test)
        END_4 = time()
        out_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('LogisticRegression', _i, train_score, test_score, ac_score, pr_score, rc_score, f_score, roc_auc_score, END_4-END_3))
        sys.stderr.write('Model: multinomial_BAYES\n')
        train_score, test_score, ac_score, pr_score, rc_score, f_score, roc_auc_score = multinomial_BAYES_func(X_train, X_test, y_train, y_test)
        END_5 = time()
        out_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('multinomial_BAYES', _i, train_score, test_score, ac_score, pr_score, rc_score, f_score, roc_auc_score, END_5-END_4))
        sys.stderr.write('Model: bernoulli_BAYES\n')
        train_score, test_score, ac_score, pr_score, rc_score, f_score, roc_auc_score = bernoulli_BAYES_func(X_train, X_test, y_train, y_test)
        END_6 = time()
        out_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('bernoulli_BAYES', _i, train_score, test_score, ac_score, pr_score, rc_score, f_score, roc_auc_score, END_6-END_5))
        sys.stderr.write('Model: gaussian_BAYES\n')
        train_score, test_score, ac_score, pr_score, rc_score, f_score, roc_auc_score = gaussian_BAYES_func(X_train, X_test, y_train, y_test)
        END_7 = time()
        out_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('gaussian_BAYES', _i, train_score, test_score, ac_score, pr_score, rc_score, f_score, roc_auc_score, END_7-END_6))
        sys.stderr.write('Model: DecisionTreeClassifier\n')
        train_score, test_score, ac_score, pr_score, rc_score, f_score, roc_auc_score = DecisionTreeClassifier_func(X_train, X_test, y_train, y_test)
        END_8 = time()
        out_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('DecisionTreeClassifier', _i, train_score, test_score, ac_score, pr_score, rc_score, f_score, roc_auc_score, END_8-END_7))
        sys.stderr.write('Model: SVM\n')
        train_score, test_score, ac_score, pr_score, rc_score, f_score, roc_auc_score = SVM_func(X_train, X_test, y_train, y_test)
        END_9 = time()
        out_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('SVM', _i, train_score, test_score, ac_score, pr_score, rc_score, f_score, roc_auc_score, END_9-END_8))
        # sys.stderr.write('Model: NN_MLPClassifier\n')
        # train_score, test_score, ac_score, pr_score, rc_score, f_score, roc_auc_score = NN_MLPClassifier_func(X_train, X_test, y_train, y_test)
        # END_10 = time()
        # out_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('NN_MLPClassifier', _i, train_score, test_score, ac_score, pr_score, rc_score, f_score, roc_auc_score, END_10-END_9))
        sys.stderr.write('Model: KNeighborsClassifier\n')
        train_score, test_score, ac_score, pr_score, rc_score, f_score, roc_auc_score = KNeighborsClassifier_func(X_train, X_test, y_train, y_test)
        END_11 = time()
        # out_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('KNeighborsClassifier', _i, train_score, test_score, ac_score, pr_score, rc_score, f_score, roc_auc_score, END_11-END_10))
        out_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('KNeighborsClassifier', _i, train_score, test_score, ac_score, pr_score, rc_score, f_score, roc_auc_score, END_11-END_9))
        # sys.exit()
    out_file.close()

def modelConstructAndPredict(inputFile, outPrefix, bcfout=None, ks_list=None, snp_list=None, test_size=0.3, model2Use=None, filter_features=False, set_divide='7,2,1', selected_features_str=None):

    if bcfout:
        vdb_median = np.median([float(_v[1]) for _k,_v in bcfout.items() if int(_v[-1]) >= 3])
        bcf_true_site = [_k for _k,_v in bcfout.items() if int(_v[-1]) >= 3 and float(_v[1]) >= vdb_median]

    # data preparation
    # inputFile = '/public/home/Songlab/xial/Work/Nanopore_RNA_editing/data/own_data/N2324962_80-1411948157_2024-02-01/results/HEK293T/own_data_HEK293T_nano_150_7/04_callSite/own_data_HEK293T_nano_a2T3_2.all.ALU.candsite.txt'
    data = pd.read_csv(inputFile, sep='\t', header=0, na_values=['nan', 'na'])
    data = data[data['REF'].isin(['A','C','G','T','a','c','g','t'])]
    data['CHROM'] = data['CHROM'].astype('str'); data['POS'] = data['POS'].astype('str')
    data['ID'] = data['CHROM'] + '_' + data['POS']
    data['MUT'] = data['REF'].str.upper()+'_'+data['ALT'].str.upper()
    # data['KNOWN'] = 2
    data.loc[data['DEL_NUM_1'] == 0, 'DEL_LEN_1'] = 0
    data.loc[data['DEL_SITE_NUM_1'] == 0, 'DEL_SITE_RATIO_1'] = 0
    data.loc[data['SKIP_1'] == 1, 'DEL_SITE_RATIO_1'] = 0
    data.loc[data['INS_NUM_1'] == 0, 'INS_LEN_1'] = 0
    data.loc[data['DEL_NUM_mid'] == 0, 'DEL_LEN_mid'] = 0
    data.loc[data['DEL_SITE_NUM_mid'] == 0, 'DEL_SITE_RATIO_mid'] = 0
    data.loc[data['SKIP_mid'] == 1, 'DEL_SITE_RATIO_mid'] = 0
    data.loc[data['INS_NUM_mid'] == 0, 'INS_LEN_mid'] = 0
    data.loc[data['DEL_NUM_3'] == 0, 'DEL_LEN_3'] = 0
    data.loc[data['DEL_SITE_NUM_3'] == 0, 'DEL_SITE_RATIO_3'] = 0
    data.loc[data['SKIP_3'] == 1, 'DEL_SITE_RATIO_3'] = 0
    data.loc[data['INS_NUM_3'] == 0, 'INS_LEN_3'] = 0
    data['SPLICE_1'] = data['SPLICE_1'] + 20
    data['SPLICE_mid'] = data['SPLICE_mid'] + 20
    data['SPLICE_3'] = data['SPLICE_3'] + 20
    data['BASE_1'] = data['MOTIF'].str[0]; data['BASE_3'] = data['MOTIF'].str[2]
    data['BASE_1_int'] = 0; data['BASE_3_int'] = 0
    data.loc[data['BASE_1'] == 'A', 'BASE_1_int'] = 1; data.loc[data['BASE_3'] == 'A', 'BASE_3_int'] = 1
    data.loc[data['BASE_1'] == 'C', 'BASE_1_int'] = 2; data.loc[data['BASE_3'] == 'C', 'BASE_3_int'] = 2
    data.loc[data['BASE_1'] == 'G', 'BASE_1_int'] = 3; data.loc[data['BASE_3'] == 'G', 'BASE_3_int'] = 3
    data.loc[data['BASE_1'] == 'T', 'BASE_1_int'] = 4; data.loc[data['BASE_3'] == 'T', 'BASE_3_int'] = 4
    data = data.drop(columns=['SEQ_ERROR_P_1', 'SEQ_ERROR_P_mid', 'SEQ_ERROR_P_3'])
    # data_snp.to_csv(outPrefix+'.snp.txt', sep='\t', index=0)
    data = data.loc[~data['ID'].isin(snp_list)] # exclude those annotated as snps
    # data_test = data[data['CHROM'].isin(['chr3', 'chr11', 'chr2', 'chr6', 'chr19', 'chrX'])]
    # data_train = data[~data['CHROM'].isin(['chr3', 'chr11', 'chr2', 'chr6', 'chr19', 'chrX'])]

    if not bcfout:
        data_known = data[data['ID'].isin(ks_list)] # all known sites
        data_known['KNOWN'] = 0
        # data_model_test.loc[(data_model_test['MUT'].isin(['A_G','T_C'])) & (data_model_test['MOTIF'].isin(['TAG', 'ATC'])), 'KNOWN'] = 1 # A>G/T>C with specific MOTIF as positive data; others as negative data
        data_known.loc[data_known['MUT'].isin(['A_G','T_C']), 'KNOWN'] = 1
        # data_known[data_known['KNOWN'] == 1].to_csv(outPrefix+'.knownSites.txt', sep='\t', index=0)
        sampled_num = min([data_known[data_known['KNOWN'] == 0].shape[0], data_known[data_known['KNOWN'] == 1].shape[0]])
        sampled_true_data = data_known[data_known['KNOWN'] == 1].sample(n=sampled_num)
        sampled_false_data = data_known[data_known['KNOWN'] == 0].sample(n=sampled_num)
    else:
        data_known = data[data['ID'].isin(list(bcfout.keys()))] # & data['ID'.isin(ks_list)]] # all known sites 
        data_known.loc[:,'KNOWN'] = 0
        data_known.loc[data_known['ID'].isin(bcf_true_site), 'KNOWN'] = 1
        sampled_num = min([data_known[data_known['KNOWN'] == 0].shape[0], data_known[data_known['KNOWN'] == 1].shape[0]])
        sampled_true_data = data_known[data_known['KNOWN'] == 1].sample(n=sampled_num)
        sampled_false_data = data_known[data_known['KNOWN'] == 0].sample(n=sampled_num)
    sys.stderr.write('#All data for training, validation and testing: %s\ntrue_num: %s; false_num: %s\n' % (data_known.shape[0], data_known[data_known['KNOWN'] == 1].shape[0], data_known[data_known['KNOWN'] == 0].shape[0]))
    data_model_test = pd.concat([sampled_true_data, sampled_false_data], axis=0)
    feature_list = list(data_model_test.iloc[:,5:57].columns.values)
    feature_list.append('BASE_1_int'); feature_list.append('BASE_3_int')
    sys.stderr.write('Sampled true_num: %s; sampled false_num: %s; number for model test: %s\n' % (sampled_true_data.shape[0], sampled_false_data.shape[0], data_model_test.shape[0]))
    data_unknown = data[~data['ID'].isin(ks_list)]
    unknown_site_list = list(data_unknown['ID'])
    sys.stderr.write('Data preparation finished.\n')

    if model2Use.split('.')[-1] != 'pkl': # to train new model
        # devide sites into training set, validation set and test set
        train_v, valid_v, test_v = [int(x) for x in set_divide.split(',')]
        train_r, valid_r, test_r = (train_v/sum([train_v, valid_v, test_v]), valid_v/sum([train_v, valid_v, test_v]), test_v/sum([train_v, valid_v, test_v]))
        data_train = data_model_test.sample(frac=train_r+valid_r)
        data_test = data_model_test.sample(frac=test_r)    
        sys.stderr.write('#data_train and #data_test: %s, %s\n' % (data_train.shape[0], data_test.shape[0]))
        sys.stderr.write('Features for analysing: '+str(feature_list)+'\n')

        X_train, X_valid, y_train, y_valid = train_test_split(data_train[feature_list], data_train['KNOWN'], test_size=test_size, random_state=0)
        if model2Use == 'RandomForestClassifier':
            Test_model, selected_features = RandomForestClassifier_func(X_train, X_valid, y_train, y_valid, outPrefix, model_selection=False, feature_list=feature_list, filter_features=filter_features, test_data_x=data_test[feature_list], test_data_y=data_test['KNOWN'])
        elif model2Use == 'HistGradientBoostingClassifier':
            Test_model, selected_features = HistGradientBoostingClassifier_func(X_train, X_valid, y_train, y_valid, outPrefix, model_selection=False, feature_list=feature_list, test_data_x=data_test[feature_list], test_data_y=data_test['KNOWN'])
        elif model2Use == 'LogisticRegression':
            Test_model, selected_features = LogisticRegression_func(X_train, X_valid, y_train, y_valid, outPrefix, model_selection=False, feature_list=feature_list, test_data_x=data_test[feature_list], test_data_y=data_test['KNOWN'])

        joblib.dump(Test_model, outPrefix+'.pkl') # store the trained model

    else:
        Test_model = joblib.load(model2Use) # to load the trained model
        if selected_features_str: 
            selected_features = selected_features_str.split(',')
            sys.stderr.write('Selected features: '+str(selected_features)+'\n')
        else:
            selected_features = feature_list
            sys.stderr.write('All features were selected.\n')
    data_unknown = data_unknown[selected_features]

    # predict unknown sites with trained model
    predict_output_df = predictUnknownSites(Test_model, data_unknown, unknown_site_list)
    predict_output_df = predict_output_df[['ID', 'Predicted_out', 'Possibility']]
    predict_output_df.to_csv(outPrefix+'.outPredict.txt', sep='\t', index=0)

    # annotate results
    # annotateResult(outPrefix+'.outPredict.txt', inputFile, ref_in, outPrefix, selected_features, annotation=None, seq_error_p=0.05)

def HistGradientBoostingClassifier_func(X_train, X_test, y_train, y_test, outPrefix, model_selection=True, feature_list=None, test_data_x=None, test_data_y=None):

    if model_selection: # to test basic scores for RandomForestClassifier, used for model comparison ['--model_selection' FUNCTION]

        test_model =  HistGradientBoostingClassifier()
        test_model.fit(X_train, y_train)
        train_score = test_model.score(X_train, y_train)
        test_score = test_model.score(X_test, y_test)
        y_pred = test_model.predict(X_test)
        ac_score, pr_score, rc_score, f_score, roc_auc_score = score(test_model, X_test, y_test)
        return(train_score, test_score, ac_score, pr_score, rc_score, f_score, roc_auc_score)
    
    else: # to train the model and predict ['--model' FUNCTION]

        from sklearn.ensemble import HistGradientBoostingClassifier
        Test_model = HistGradientBoostingClassifier()
        Test_model.fit(X_train, y_train)
        sys.stderr.write("Results of RF fitting: \n")
        sys.stderr.write('Accuracy on training set: %s\n' % (Test_model.score(X_train, y_train)))
        sys.stderr.write('Accuracy on test set: %s\n' % (Test_model.score(X_test, y_test)))

        # model evaluation 
        sys.stderr.write('HGB Training set:\n')
        predict_Target=Test_model.predict(X_train)
        sys.stderr.write(str(metrics.classification_report(y_train,predict_Target))+'\n')
        sys.stderr.write(str(metrics.confusion_matrix(y_train, predict_Target))+'\n')
        sys.stderr.write('HGB Validation set:\n')
        predict_target=Test_model.predict(X_test)
        sys.stderr.write(str(metrics.classification_report(y_test,predict_target))+'\n')
        sys.stderr.write(str(metrics.confusion_matrix(y_test, predict_target))+'\n')

        predict_target_prob = Test_model.predict_proba(test_data_x) # calculate the possibility for each cluster (0, 1)
        predict_target_prob_rf = predict_target_prob[:,1] # the possibility of predicted as "1"
        predict_rf = Test_model.predict(test_data_x) # to output prediction
        predict_output_df = pd.DataFrame({'prob_1':predict_target_prob_rf,'predict_out': predict_rf,'true_out':list(test_data_y)})
        predict_output_df.to_csv(outPrefix+'.test_set_predict.txt', sep='\t', index=0)
        ROC_plot(list(test_data_y), list(predict_target_prob_rf), outPrefix)
        
        return(Test_model, feature_list) # Output the trained model and the selected features

def RandomForestClassifier_func(X_train, X_test, y_train, y_test, outPrefix, model_selection=True, feature_list=None, filter_features=False, test_data_x=None, test_data_y=None):

    if model_selection: # to test basic scores for RandomForestClassifier, used for model comparison ['--model_selection' FUNCTION]

        test_model =  RandomForestClassifier()
        test_model.fit(X_train, y_train)
        train_score = test_model.score(X_train, y_train)
        test_score = test_model.score(X_test, y_test)
        y_pred = test_model.predict(X_test)
        ac_score, pr_score, rc_score, f_score, roc_auc_score = score(test_model, X_test, y_test)
        return(train_score, test_score, ac_score, pr_score, rc_score, f_score, roc_auc_score)
    
    else: # to train the model and predict ['--model' FUNCTION]

        from sklearn.ensemble import RandomForestClassifier
        test_model = RandomForestClassifier()
        test_model.fit(X_train, y_train)
        sys.stderr.write("Results of RF fitting: \n")
        sys.stderr.write('Accuracy on training set: %s\n' % (test_model.score(X_train, y_train)))
        sys.stderr.write('Accuracy on test set: %s\n' % (test_model.score(X_test, y_test)))

        if filter_features:
            # feature importance evaluation
            selected_features = featureSelection(X_train, test_model, feature_list, outPrefix)
            X_train = X_train[selected_features]
            sys.stderr.write('Selected features: '+','.join(list(X_train.columns.values))+'\n')
            test_model_selected = RandomForestClassifier() # use selected features to train new model in case overfitting
            test_model_selected.fit(X_train, y_train)
            sys.stderr.write("Results of RF fitting after feature selection: \n")
            sys.stderr.write('Accuracy on training set: %s\n' % (test_model_selected.score(X_train, y_train)))
            sys.stderr.write('Accuracy on test set: %s\n' % (test_model_selected.score(X_test[selected_features], y_test)))
            Test_model = test_model_selected
            X_test = X_test[selected_features]; test_data_x = test_data_x[selected_features] # feature filtration
        else:
            Test_model = test_model; selected_features = feature_list.copy()

        # model evaluation 
        sys.stderr.write('RF Training set:\n')
        predict_Target=Test_model.predict(X_train)
        sys.stderr.write(str(metrics.classification_report(y_train,predict_Target))+'\n')
        sys.stderr.write(str(metrics.confusion_matrix(y_train, predict_Target))+'\n')
        sys.stderr.write('RF Validation set:\n')
        predict_target=Test_model.predict(X_test)
        sys.stderr.write(str(metrics.classification_report(y_test,predict_target))+'\n')
        sys.stderr.write(str(metrics.confusion_matrix(y_test, predict_target))+'\n')

        predict_target_prob = Test_model.predict_proba(test_data_x) # calculate the possibility for each cluster (0, 1)
        predict_target_prob_rf = predict_target_prob[:,1] # the possibility of predicted as "1"
        predict_rf = Test_model.predict(test_data_x) # to output prediction
        predict_output_df = pd.DataFrame({'prob_1':predict_target_prob_rf,'predict_out': predict_rf,'true_out':list(test_data_y)})
        predict_output_df.to_csv(outPrefix+'.test_set_predict.txt', sep='\t', index=0)
        ROC_plot(list(test_data_y), list(predict_target_prob_rf), outPrefix)
        
        return(Test_model, selected_features) # Output the trained model and the selected features

def predictUnknownSites(Test_model, unknown_data_x, site_list):

    predict_target_prob = Test_model.predict_proba(unknown_data_x) # calculate the possibility for each cluster (0, 1)
    predict_target_prob_rf = predict_target_prob[:,1] # the possibility of predicted as "1"
    predict_rf = Test_model.predict(unknown_data_x) # to output prediction
    predict_output_df = pd.DataFrame({'Predicted_out': predict_rf, 'Possibility':predict_target_prob_rf, 'ID':site_list})
    return(predict_output_df)

def featureSelection(X_train, test_model, feature_list, outPrefix):

    # feature importance evaluation
    importances = test_model.feature_importances_
    indices = np.argsort(importances)[::-1]
    sys.stderr.write('Feature_importances_: \n')
    importances_dict = {}
    for f in range(X_train.shape[1]):
        importances_dict.setdefault('weight', []).append(importances[indices[f]])
        importances_dict.setdefault('Feature', []).append(feature_list[indices[f]])
        sys.stderr.write("%2d) %-*s\t%f\n" % (f + 1, 30, feature_list[indices[f]], importances[indices[f]])) # Output feature importances
    importances_df = pd.DataFrame(importances_dict)
    importances_df.to_csv(outPrefix+'.feature_importance.txt', sep='\t', index=0)
    featureImportance_plot(importances_df, outPrefix)
    importances_value = list(importances_df['weight'])
    # thresh = np.median(importances_value[:int(len(importances_value)*0.5)]) # Third quartile for importances
    thresh = np.median(importances_value)
    selected_features = list(importances_df.loc[importances_df['weight'] > thresh, 'Feature']) # selected features

    return(selected_features)

def LinearRegression_func(X_train, X_test, y_train, y_test):

    test_model =  LinearRegression()
    test_model.fit(X_train, y_train)
    train_score = test_model.score(X_train, y_train)
    test_score = test_model.score(X_test, y_test)
    y_pred = test_model.predict(X_test)
    ac_score, pr_score, rc_score, f_score, roc_auc_score = score(test_model, X_test, y_test)
    return(train_score, test_score, ac_score, pr_score, rc_score, f_score, roc_auc_score)

def LogisticRegression_func(X_train, X_test, y_train, y_test, outPrefix, model_selection=True, feature_list=None, test_data_x=None, test_data_y=None):

    if model_selection:

        test_model = LogisticRegression()
        test_model.fit(X_train, y_train)
        train_score = test_model.score(X_train, y_train)
        test_score = test_model.score(X_test, y_test)
        y_pred = test_model.predict(X_test)
        ac_score, pr_score, rc_score, f_score, roc_auc_score = score(test_model, X_test, y_test)
        return(train_score, test_score, ac_score, pr_score, rc_score, f_score, roc_auc_score)

    else:

        from sklearn.linear_model import LogisticRegression
        Test_model = LogisticRegression()
        Test_model.fit(X_train, y_train)
        sys.stderr.write("Results of RF fitting: \n")
        sys.stderr.write('Accuracy on training set: %s\n' % (Test_model.score(X_train, y_train)))
        sys.stderr.write('Accuracy on test set: %s\n' % (Test_model.score(X_test, y_test)))

        # model evaluation 
        sys.stderr.write('LogisticR Training set:\n')
        predict_Target=Test_model.predict(X_train)
        sys.stderr.write(str(metrics.classification_report(y_train,predict_Target))+'\n')
        sys.stderr.write(str(metrics.confusion_matrix(y_train, predict_Target))+'\n')
        sys.stderr.write('LogisticR Validation set:\n')
        predict_target=Test_model.predict(X_test)
        sys.stderr.write(str(metrics.classification_report(y_test,predict_target))+'\n')
        sys.stderr.write(str(metrics.confusion_matrix(y_test, predict_target))+'\n')

        predict_target_prob = Test_model.predict_proba(test_data_x) # calculate the possibility for each cluster (0, 1)
        predict_target_prob_rf = predict_target_prob[:,1] # the possibility of predicted as "1"
        predict_rf = Test_model.predict(test_data_x) # to output prediction
        predict_output_df = pd.DataFrame({'prob_1':predict_target_prob_rf,'predict_out': predict_rf,'true_out':list(test_data_y)})
        predict_output_df.to_csv(outPrefix+'.test_set_predict.txt', sep='\t', index=0)
        ROC_plot(list(test_data_y), list(predict_target_prob_rf), outPrefix)
        
        return(Test_model, feature_list) # Output the trained model and the selected features

def multinomial_BAYES_func(X_train, X_test, y_train, y_test):

    #test_model = sk_bayes.MultinomialNB(alpha=1.0,fit_prior=True,class_prior=None)
    test_model = sk_bayes.MultinomialNB()
    test_model.fit(X_train, y_train)
    train_score = test_model.score(X_train, y_train)
    test_score = test_model.score(X_test, y_test)
    y_pred = test_model.predict(X_test)
    ac_score, pr_score, rc_score, f_score, roc_auc_score = score(test_model, X_test, y_test)
    return(train_score, test_score, ac_score, pr_score, rc_score, f_score, roc_auc_score)

def bernoulli_BAYES_func(X_train, X_test, y_train, y_test):

    #test_model = sk_bayes.BernoulliNB(alpha=1.0,binarize=0.0,fit_prior=True,class_prior=None)
    test_model = sk_bayes.BernoulliNB()
    test_model.fit(X_train, y_train)
    train_score = test_model.score(X_train, y_train)
    test_score = test_model.score(X_test, y_test)
    y_pred = test_model.predict(X_test)
    ac_score, pr_score, rc_score, f_score, roc_auc_score = score(test_model, X_test, y_test)
    return(train_score, test_score, ac_score, pr_score, rc_score, f_score, roc_auc_score)

def gaussian_BAYES_func(X_train, X_test, y_train, y_test):

    test_model = sk_bayes.GaussianNB()
    test_model.fit(X_train, y_train)
    train_score = test_model.score(X_train, y_train)
    test_score = test_model.score(X_test, y_test)
    y_pred = test_model.predict(X_test)
    ac_score, pr_score, rc_score, f_score, roc_auc_score = score(test_model, X_test, y_test)
    return(train_score, test_score, ac_score, pr_score, rc_score, f_score, roc_auc_score)

def DecisionTreeClassifier_func(X_train, X_test, y_train, y_test):

    #test_model = DecisionTreeClassifier(criterion='entropy',max_depth=None,min_samples_split=2,min_samples_leaf=1,max_features=None,max_leaf_nodes=None,min_impurity_decrease=0)
    test_model = DecisionTreeClassifier()
    test_model.fit(X_train, y_train)
    train_score = test_model.score(X_train, y_train)
    test_score = test_model.score(X_test, y_test)
    y_pred = test_model.predict(X_test)
    ac_score, pr_score, rc_score, f_score, roc_auc_score = score(test_model, X_test, y_test)
    return(train_score, test_score, ac_score, pr_score, rc_score, f_score, roc_auc_score)

def SVM_func(X_train, X_test, y_train, y_test):

    #test_model = SVC(C=1.0,kernel='rbf',gamma='auto')
    test_model = SVC()
    test_model.fit(X_train, y_train)
    train_score = test_model.score(X_train, y_train)
    test_score = test_model.score(X_test, y_test)
    y_pred = test_model.predict(X_test)
    ac_score, pr_score, rc_score, f_score, roc_auc_score = score(test_model, X_test, y_test)
    return(train_score, test_score, ac_score, pr_score, rc_score, f_score, roc_auc_score)

def NN_MLPClassifier_func(X_train, X_test, y_train, y_test):

    #test_model = MLPClassifier(activation='tanh',solver='adam',alpha=0.0001,learning_rate='adaptive',learning_rate_init=0.001,max_iter=200)
    test_model = MLPClassifier()
    test_model.fit(X_train, y_train)
    train_score = test_model.score(X_train, y_train)
    test_score = test_model.score(X_test, y_test)
    y_pred = test_model.predict(X_test)
    ac_score, pr_score, rc_score, f_score, roc_auc_score = score(test_model, X_test, y_test)
    return(train_score, test_score, ac_score, pr_score, rc_score, f_score, roc_auc_score)

def KNeighborsClassifier_func(X_train, X_test, y_train, y_test):

    test_model = KNeighborsClassifier() 
    test_model.fit(X_train, y_train)
    train_score = test_model.score(X_train.values, y_train.values)
    test_score = test_model.score(X_test.values, y_test.values)
    # sys.stderr.write('train_score: %s; test_score: %s\n' % (train_score, test_score))
    y_pred = test_model.predict(X_test.values)
    # sys.stderr.write('y_pred: '+str(y_pred[:10])+'\n')
    ac_score, pr_score, rc_score, f_score, roc_auc_score = score(test_model, X_test, y_test, model='KNeighborsClassifier')
    return(train_score, test_score, ac_score, pr_score, rc_score, f_score, roc_auc_score)

def score(test_model, X_test, y_test, model=None):

    if model != 'KNeighborsClassifier':
        ac_score = np.mean(cross_val_score(test_model, X_test, y=y_test, scoring='accuracy', cv=10, n_jobs=1))
        pr_score = np.mean(cross_val_score(test_model, X_test, y=y_test, scoring='precision', cv=10, n_jobs=1))
        rc_score = np.mean(cross_val_score(test_model, X_test, y=y_test, scoring='recall', cv=10, n_jobs=1))
        f_score = np.mean(cross_val_score(test_model, X_test, y=y_test, scoring='f1', cv=10, n_jobs=1))
        roc_auc_score = np.mean(cross_val_score(test_model, X_test, y=y_test, scoring='roc_auc', cv=10, n_jobs=1))
    else:
        ac_score = np.mean(cross_val_score(test_model, X_test.values, y=y_test.values, scoring='accuracy', cv=10, n_jobs=1))
        pr_score = np.mean(cross_val_score(test_model, X_test.values, y=y_test.values, scoring='precision', cv=10, n_jobs=1))
        rc_score = np.mean(cross_val_score(test_model, X_test.values, y=y_test.values, scoring='recall', cv=10, n_jobs=1))
        f_score = np.mean(cross_val_score(test_model, X_test.values, y=y_test.values, scoring='f1', cv=10, n_jobs=1))
        roc_auc_score = np.mean(cross_val_score(test_model, X_test.values, y=y_test.values, scoring='roc_auc', cv=10, n_jobs=1))
    return(ac_score, pr_score, rc_score, f_score, roc_auc_score)

def featureImportance_plot(importances_df, outPrefix):

    df=importances_df.sort_values(by='weight')
    data_height=df['weight'].values.tolist()
    data_x=df['Feature'].values.tolist()

    font = {'family': 'ARIAL','size': 7}
    # sns.set(font_scale=1.2)
    plt.rc('font',family='ARIAL')

    plt.figure(figsize=(8,10))
    plt.barh(range(len(data_x)), data_height, color='#6699CC')
    plt.yticks(range(len(data_x)),data_x,fontsize=8)

    plt.tick_params(labelsize=8)
    plt.xlabel('Feature importance',fontsize=10)
    plt.title("Feature importance analysis",fontsize = 10)
    plt.savefig(outPrefix+'.feature_imp.pdf', )
    # plt.show()

def ROC_plot(list1, list2, outPrefix):

    fpr_model, tpr_model, thresholds = roc_curve(list1, list2, pos_label=1)
    roc_auc_model = auc(fpr_model,tpr_model)

    font = {'family': 'ARIAL','size': 12}
    # sns.set(font_scale=1.2)
    plt.rc('font',family='ARIAL')
    plt.figure(figsize=(8,9))

    plt.plot(fpr_model,tpr_model,'blue',label='AUC = %0.2f'% roc_auc_model)
    plt.legend(loc='lower right',fontsize = 12)
    plt.plot([0,1],[0,1],'r--')
    plt.title("ROC curve for test set (No.=%s)" % (len(list1)), fontsize = 14)
    plt.ylabel('True Positive Rate',fontsize = 14)
    plt.xlabel('Flase Positive Rate',fontsize = 14)
    plt.savefig(outPrefix+'.ROC_curve.pdf')
    # plt.show()

def annotateResult(inputFile, cand_site_info, ref_in, out_prefix, selected_features, annotation=None, seq_error_p=0.05):

    true_anno_out_file = open(out_prefix+'.trueAnno.txt', 'w')
    true_anno_out_file.write('ID\tREF\tALT\tREADS_NUM\tEDIT_LEVEL\tMOTIF_SEQ\tPROB\tSEQ_ERROR_P\t'+'\t'.join(selected_features)+'\tGENE_ID\tTRANS_ID\n'+'\n')

    a = 0; mut_list = []
    for line in open(inputFile, 'r'):
        if line.startswith('ID'): continue
        a += 1
        if a % 10000 == 0: sys.stderr.write('%s candidate sites were analyzed.\n' % (a))
        info = line.strip().split('\t')
        pos_id = info[0]; true_false = info[1]; score = info[2]
        # if pos_id != 'chr1_10519167': continue
        if true_false != '1': continue # only output predicted TRUE sites
        
        chrom = pos_id.split('_')[0]; pos = int(pos_id.split('_')[1])
        trans_id = []; gene_id = []; strand_list = []
        for _i in annotation[chrom]:
            if _i[1][0] < pos <= _i[1][1]: # pos in the target transcript
                strand_list.append(_i[4]); gene_id.append(_i[5][0]); trans_id.append(_i[0])
            elif pos < _i[1][0]: break # no more possible genes
            else: continue
        if len(set(strand_list)) == 1: strand = strand_list[0]
        else: strand = None
        if trans_id: trans_id_list = ','.join(trans_id)
        else: trans_id_list = 'None'
        if gene_id: gene_id_list = ','.join(gene_id)
        else: gene_id_list = 'None'

        if not strand: ID = '%s_%s_%s' % (chrom, None, pos)
        else: ID = '%s_%s_%s' % (chrom, strand, pos)
        if strand == '+': motif_seq = ref_in.fetch(chrom, pos-3, pos+2).upper()
        elif strand == '-': motif_seq = ''.join([base_pair(x) for x in ref_in.fetch(chrom, pos-3, pos+2)][::-1]).upper()
        else: motif_seq = ref_in.fetch(chrom, pos-3, pos+2).upper()
        _pos, reads_num, ref, alt, edit_level, depth = [x for x in cand_site_info[chrom] if x[0] == pos][0]
        seq_error_P = round(binom.pmf(k=int(reads_num.split(',')[1]), n=depth, p=seq_error_p), 4)
        mut_list.append(ref.upper()+'_'+alt.upper())
        true_anno_out_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (ID, ref.upper(), alt.upper(), reads_num, edit_level, motif_seq, score, seq_error_P, gene_id_list, trans_id_list))
    
    true_anno_out_file.close()

    # output the mut_ratio
    mut_count = Counter(mut_list)
    for _m, _c in sorted(mut_count.items(), key=lambda x:x[0]):
        sys.stderr.write('%s\t%s\t%s\n' % (_m, _c, _c/sum(mut_count.values())))

if __name__ == '__main__':
    arg(sys.argv)
