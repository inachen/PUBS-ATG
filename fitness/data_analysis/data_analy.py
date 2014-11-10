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

UBIQ_SEQ = "MQIFVKTLTGKTITLEVESSDTIDNVKSKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG!"
AA_LST= ["STOP", "W", "F", "Y", "L", "I", "M", "V", "C", "A", "G", "P", "S", "T", "N", "Q", "H", "R", "K", "D", "E"]

# advantageous mutations



def add_timecourses(t1, t2):
    t1d1 = t1[0]
    t1d2 = t1[1]
    t2d1 = t2[0]
    t2d2 = t2[1]
    t = ([t1d1[0] + t2d1[0], t1d1[1] + t2d1[1], t1d1[2] + t2d1[2]], 
        [t1d2[0] + t2d2[0], t1d2[1] + t2d2[1], t1d2[2] + t2d2[2]])

    return t

def ratio_slope(lst1, lst2):

    # lst1/lst2
    p1 = lst1[0]/float(lst2[0])
    p2 = lst1[2]/float(lst2[2])

    # return t2 and t3 sel coef
    # s1 = lst[1]-lst[]
    return ()



def calc_sel_coeff():
    ### INCOMPLETE FUNCTION ###

    # open time course data
    infile = open("input_data/filtered_seq_perfect.pkl", 'rb')
    time_dic = pickle.load(infile)
    infile.close()

    # wt keys
    wt_seq = [(i, UBIQ_SEQ[i-1]) for i in range (1,77)]

    # holder for wt counts
    wt_timecourse = ([0,0,0], [0,0,0])

    mut_time_dic = copy.deepcopy(time_dic)

    # calculate total WT counts
    for mut, times in time_dic.iteritems():

        # check if mutation is wt
        if mut in wt_seq:
            print mut
            mut_timecourse = add_timecourses(times, wt_timecourse)
            mut_time_dic.pop(mut, None)

        elif mut == (0, 'WT'):
            print "wt"
            mut_timecourse = add_timecourses(times, wt_timecourse)
            mut_time_dic.pop(mut, None)

    # record remaining mutations
    all_muts = mut_time_dic.keys()

    sel_coef_mat = np.zeros((21,77))

    # use -1 for wt
    sel_coef_mat.fill(-1)

    print sel_coef_mat

    # calulate selection coefficients
    # for mut, times in mut_time_dic.iteritems():

        # calculate ratio slope (mut/wt)
        # r1s = []


    # export to csv
    # finaldata = pd.DataFrame(, index=names, columns=names)
    # df.to_csv('df.csv', index=True, header=True, sep=' ')

def get_consensus():
    # get all the amino acid sequences

    # open file
    infile = open("input_data/filtered_seq_perfect.pkl", 'rb')
    t_dic = pickle.load(infile)
    infile.close()

    ubiq_lst = list(UBIQ_SEQ)

    # track all seqs
    seq_d1_lst = []
    seq_d2_lst = []

    # wt counts
    count_1_1 = 0
    count_1_3 = 0
    count_2_1 = 0
    count_2_3 = 0

    time_dic = copy.deepcopy(t_dic)

    # wt keys
    wt_seq = [(i, UBIQ_SEQ[i-1]) for i in range (1,77)]

    # calculate total WT counts
    for mut, times in t_dic.iteritems():

        # check if mutation is wt
        if mut in wt_seq:
            time_dic.pop(mut, None)

        elif mut == (0, 'WT'):
            time_dic.pop(mut, None)

    # get wt counts at each time
    for mut, times in time_dic.iteritems():

        d1_1 = times[0][0]
        d1_3 = times[0][2]
        d2_1 = times[1][0]
        d2_3 = times[1][2]

        count_1_1 += d1_1
        count_1_3 += d1_3
        count_2_1 += d2_1
        count_2_3 += d2_3

    for mut, times in time_dic.iteritems():

        mut_d1_1 = times[0][0]
        mut_d1_3 = times[0][2]
        mut_d2_1 = times[1][0]
        mut_d2_3 = times[1][2]

        ratio_d1_1 = int(mut_d1_1/float(count_1_1) * 10000)
        ratio_d1_3 = int(mut_d1_3/float(count_1_3) * 10000)
        ratio_d2_1 = int(mut_d2_1/float(count_2_1) * 10000)
        ratio_d2_3 = int(mut_d2_3/float(count_2_3) * 10000)

        pos = mut[0]
        aa = mut[1]

        # seq = copy.deepcopy(ubiq_lst)
        seq = ['-' for i in range(77)]

        if pos != 0:
            if aa != "STOP":
                seq[pos-1] = aa
            else :
                seq[pos-1] = "!"
        
        if ratio_d1_3 > ratio_d1_1:
            seq_d1_lst.extend(["".join(seq) for i in range(ratio_d1_3)])

        if ratio_d2_3 > ratio_d2_1:
            seq_d2_lst.extend(["".join(seq) for i in range(ratio_d2_3)])

    cw = csv.writer(open("consensus_seq_d1.csv", "wb"))

    for s in seq_d1_lst:
        cw.writerow([s])

    cw2 = csv.writer(open("consensus_seq_d2.csv", "wb"))

    for s in seq_d2_lst:
        cw2.writerow([s])

get_consensus()


