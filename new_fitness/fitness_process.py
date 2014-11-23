import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import collections
from collections import namedtuple
import csv
import numpy as np
import math
import matplotlib.cm as cm
from mpltools import style
from mpltools import color
import copy
import os.path
import scipy.stats as stats

from glob import glob
from Bio import SeqIO
import cPickle as pickle

UBIQ_SEQ = "MQIFVKTLTGKTITLEVESSDTIDNVKSKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG!"
AA_LST= ["STOP", "W", "F", "Y", "L", "I", "M", "V", "C", "A", "G", "P", "S", "T", "N", "Q", "H", "R", "K", "D", "E"]
POS_NUM = 77

# arbitrary complex order
NAMES = ['DMSO', 'Caffeine', 'Hydroxyurea']

# load data pickles
dmso_in = open('data_pkls/Codon_FitScore_DMSO_day1.pkl', 'rb')
dmso_codon = pickle.load(dmso_in)
dmso_in.close()

caffeine_in = open('data_pkls/Codon_FitScore_Caffeine_day2.pkl', 'rb')
caffeine_codon = pickle.load(caffeine_in)
caffeine_in.close()

hu_in = open('data_pkls/Codon_FitScore_HU_day1.pkl', 'rb')
hu_codon = pickle.load(hu_in)
hu_in.close()

# load codon to amino acid dictionary
trans_in = open('dic_pkls/translate.pkl', 'rb')
translate = pickle.load(trans_in)
trans_in.close()

# pos, aa -> fitness dictionaries
dmso_fit = {}
caffeine_fit = {}
hu_fit = {}

codon_dic_lst = [dmso_codon, caffeine_codon, hu_codon]
fit_dic_lst = [dmso_fit, caffeine_fit, hu_fit]

for dic_codon, dic_aa in zip(codon_dic_lst, fit_dic_lst):
    for k, val in dic_codon.iteritems():

        new_k = k

        if k[1] != 'WT':
            new_k = (k[0], translate[k[1].replace('T', 'U')])

        dic_aa[new_k] = val[0]

# matrix of fitness values
dmso_fit_mat = np.zeros((len(AA_LST), POS_NUM))
dmso_fit_mat = dmso_fit_mat - 1
caffeine_fit_mat = np.zeros((len(AA_LST), POS_NUM))
caffeine_fit_mat = caffeine_fit_mat - 1
hu_fit_mat = np.zeros((len(AA_LST), POS_NUM))
hu_fit_mat = hu_fit_mat - 1

fit_mat_lst = [dmso_fit_mat, caffeine_fit_mat, hu_fit_mat]

# generate heatmaps
for num, fit_mat in enumerate(fit_mat_lst):
    for i in range(1, POS_NUM+1):
        for j, aa in enumerate(AA_LST):
            row = j
            col = i - 1

            key = (i, aa)

            if key in fit_dic_lst[num]:
                fit_mat[row, col] = fit_dic_lst[num][key]


    np.savetxt("fitness_csv/" + NAMES[num] + "_fitness.csv", fit_mat, delimiter=",")


style.use("ggplot")

# generate difference heatmaps
for n1, mat1 in enumerate(fit_mat_lst):

    for n2, mat2 in enumerate(fit_mat_lst):

        if n1 != n2:
            # difference csv
            sub_mat = np.subtract(mat1, mat2)

            np.savetxt("fitness_csv/" + NAMES[n1] + '_sub_' + NAMES[n2] + "_fitness.csv", 
                sub_mat, delimiter=",")

            # correlation plot
            mat1_lst = np.array(mat1).reshape(-1,).tolist()
            mat2_lst = np.array(mat2).reshape(-1,).tolist()

            # m, b = np.polyfit(mat1_lst, mat2_lst, 1)

            m, b, r_value, p_value, std_err = stats.linregress(mat1_lst, mat2_lst)

            r_sq = r_value**2

            mat1_lst_fit = [x * float(m) + b for x in mat1_lst]

            fig = plt.figure()
            plt.plot(mat1_lst, mat2_lst, '.', color=plt.rcParams['axes.color_cycle'][1])
            plt.plot(mat1_lst, mat1_lst_fit, color='black')
            plt.annotate("r2 = %.3f" % r_sq, xy=(0.05, 0.95), xycoords='axes fraction', fontsize=16,
                xytext=(0.5, -1), textcoords='offset points',
                ha='left', va='top')
            plt.xlabel(NAMES[n1])
            plt.ylabel(NAMES[n2])
            fig.savefig("plots/" + NAMES[n1] + '_vs_' + NAMES[n2] + "_corr.png")

# generate txt files for sequence logos

# just HU
# sequences more fit than wildtype

ubq_lst = list(UBIQ_SEQ)

HU_seq_lst = []

# for storing prop dic
HU_fitness_prop_dic = {}

HU_sub_DMSO_prop_dic = {}

for key, val in hu_fit.iteritems():

    if val > 0 and key[1] != 'WT' and key[1] != 'STOP':

        num = val * 100

        seq_lst = ['-' for i in range(len(ubq_lst))]
        seq_lst[key[0]-1] = key[1]

        for i in range(int(num)):
            HU_seq_lst.append("".join(seq_lst))


cw = csv.writer(open("seq_logo/HU_consensus.csv", "wb"))

for s in HU_seq_lst:
    cw.writerow([s])

# HU minus DMSO
HU_sub_dmso_seq_lst = []

hu_sub_dmso_mat = np.subtract(hu_fit_mat, dmso_fit_mat)

hu_sub_dmso_unfit_seq_lst = []

for (i,j), value in np.ndenumerate(hu_sub_dmso_mat):

    if j == 62:
        print AA_LST[i]
        print value

    if value > 0:

        if ubq_lst[j] != AA_LST[i]:

            if AA_LST[i] != 'STOP':

                num = value * 100

                seq_lst = ['-' for x in range(len(ubq_lst))]
                seq_lst[j] = AA_LST[i]

                for i in range(int(num)):
                    HU_sub_dmso_seq_lst.append("".join(seq_lst))
 
    elif value < 0:

            val2 = -value

            if ubq_lst[j] != AA_LST[i]:

                if AA_LST[i] != 'STOP':

                    num = val2 * 100

                    seq_lst = ['-' for x in range(len(ubq_lst))]
                    seq_lst[j] = AA_LST[i]

                    for i in range(int(num)):
                        hu_sub_dmso_unfit_seq_lst.append("".join(seq_lst))

cw = csv.writer(open("seq_logo/HU_sub_DMSO_consensus.csv", "wb"))

for s in HU_sub_dmso_seq_lst:
    cw.writerow([s])      

cw = csv.writer(open("seq_logo/hu_sub_DMSO_unfit_consensus.csv", "wb"))

for s in hu_sub_dmso_unfit_seq_lst:
    cw.writerow([s])            






            






