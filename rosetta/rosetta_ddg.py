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
# from mpltools import style
# from mpltools import color
from Bio import SeqIO
import cPickle as pickle
import pandas as pd
import json

# Constants

# file paths
IN_DIR = "input_jsons"
OUT_DIR = "out_pickles"
FSEP = "/"
FP_MONO = "uby_1ubq.json"

# parsing mutations A X1Y, 1 is position, Y is mutation
POS_POS = 3
MUT_POS = 4

# polarity binning
HYDROPHOBIC = ['M', 'C', 'I', 'L', 'Y', 'F', 'W']
POLAR = ['Q', 'P', 'N', 'A', 'T', 'S', 'V', 'G']
CHARGED = ['E', 'D', 'K', 'R', 'H']

def json_to_dic(jfile, param):

    '''convert json file of Rosetta DDG values to dictionary'''

    json_data = open(jfile, 'r')
    raw_ddg = json.load(json_data)

    ddg_dic = {}

    for dat in raw_ddg['data']:

        mut = dat['Mutation']

        # get the position and amino acid
        if len(mut) == 5:
            allel = (int(mut[POS_POS]), str(mut[MUT_POS]))
        else:
            allel = (int(mut[POS_POS:POS_POS+2]), str(mut[MUT_POS+1]))

        ddg_dic[allel] = dat[param]

    return ddg_dic

def bin_by_hydrophobe(ddg_dic):
    '''bin by hydrophobicity'''

    hydro_bin_dic = {'hydrophobic': {}, 'polar':{}, 'charged':{}}

    for k, val in ddg_dic.iteritems():
        aa = k[1]

        if aa in HYDROPHOBIC:
            hydro_bin_dic['hydrophobic'][k] = val
        elif aa in POLAR:
            hydro_bin_dic['polar'][k] = val
        else:
            hydro_bin_dic['charged'][k] = val

    return hydro_bin_dic

def dic_to_lsts(dic):
    # TODO

def run():
    # parameter value to extract
    param = "global_DDG"

    # get ddg dic
    ddg_dic = json_to_dic(IN_DIR + FSEP + FP_MONO, param)

    pickle.dump(ddg_dic, open(OUT_DIR + FSEP + 'ddg_dic.pkl', 'wb'))

    # bin by hydrophobicity
    hydro_bin_dic = bin_by_hydrophobe(ddg_dic)

    pickle.dump(ddg_dic, open(OUT_DIR + FSEP + 'hydro_bin_dic.pkl', 'wb'))


run()