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
import copy
import os.path

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
caffeine_fit_mat = np.zeros((len(AA_LST), POS_NUM))
hu_fit_mat = np.zeros((len(AA_LST), POS_NUM))

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

# generate difference heatmaps
for n1, mat1 in enumerate(fit_mat_lst):

    for n2, mat2 in enumerate(fit_mat_lst):

        if n1 != n2:
            sub_mat = np.subtract(mat1, mat2)

            np.savetxt("fitness_csv/" + NAMES[n1] + '_sub_' + NAMES[n2] + "_fitness.csv", 
                sub_mat, delimiter=",")

            






