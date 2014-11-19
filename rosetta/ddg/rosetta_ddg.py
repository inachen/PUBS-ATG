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
IN_DIR = "input"
OUT_DIR = "out_pickles"
FSEP = "/"
FP_MONO = "uby_1ubq.json"
FP_CUE = "uby_cue.json"
FP_OTU = "uby_otu.json"
FP_RPN13 = "uby_rpn13.json"
FP_SH3 = "uby_sh3.json"
FP_UQCON = "uby_uqcon.json"
FITNESS_PKL = "D2S3fitness_scores.pkl"

BINDING_LST = [("Monomer", FP_MONO), ("CUE", FP_CUE), ("OTU", FP_OTU),
    ("RPN13", FP_RPN13), ("SH3", FP_SH3), ("UQCON", FP_UQCON)]

BURIED = [3, 5, 15, 17, 23, 26, 30, 41, 43, 45, 50, 56, 59, 61, 67, 69]

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

def bin_by_topo(ddg_dic):
    topo_bin_dic = {'buried': {}, 'surface':{}}

    for k, val in ddg_dic.iteritems():
        pos = k[0]

        if pos in BURIED:
            topo_bin_dic['buried'][k] = val
        else:
            topo_bin_dic['surface'][k] = val

    return topo_bin_dic

def bin_dic_to_csv(dic, fit_dic, name):
    
    for label, d in dic.iteritems():
        ddg_lst, fit_lst = dic_to_lsts(d, fit_dic)

        lsts_to_csv(ddg_lst, fit_lst, 'DDG', 'Fitness', label, name)

def dic_to_lsts(dic, fit_dic):

    ddg_lst = []
    fit_lst = []

    for mut, ddg in dic.iteritems():

        if mut in fit_dic.keys():
            ddg_lst.append(ddg)
            fit_lst.append(fit_dic[mut])

    return (ddg_lst, fit_lst)

def lsts_to_csv(lst1, lst2, h1, h2, bintype, name):

    lst = [list(i) for i in zip(lst1, lst2)]

    outfile = open(OUT_DIR + FSEP + bintype+"_"+name+"_"+h1+"_"+h2+".csv",'wb')

    wr = csv.writer(outfile)
    wr.writerow([h1, h2])
    for i in lst:
        wr.writerow(i)

def run():
    # parameter value to extract
    param = "global_DDG"

    # get fitness dictionary
    fitness_dic = pickle.load(open(IN_DIR+FSEP+FITNESS_PKL, "rb"))

    for (name, fp) in BINDING_LST:

        # get ddg dic
        ddg_dic = json_to_dic(IN_DIR+FSEP+fp, param)

        pickle.dump(ddg_dic, open(OUT_DIR+FSEP+'ddg_dic_'+name+'.pkl', 'wb'))

        # bin by hydrophobicity
        hydro_bin_dic = bin_by_hydrophobe(ddg_dic)

        pickle.dump(ddg_dic, open(OUT_DIR+FSEP+'hydro_bin_dic_'+name+'.pkl', 'wb'))

        # bin by topology
        topo_bin_dic = bin_by_topo(ddg_dic)

        pickle.dump(ddg_dic, open(OUT_DIR+FSEP+'topo_bin_dic_'+name+'.pkl', 'wb'))

        # generate lists for corrleation
        ddg_lst, fit_lst = dic_to_lsts(ddg_dic, fitness_dic)

        lsts_to_csv(ddg_lst, fit_lst, 'DDG', 'Fitness', 'all', name)

        bin_dic_to_csv(hydro_bin_dic, fitness_dic, name)
        bin_dic_to_csv(topo_bin_dic, fitness_dic, name)


run()