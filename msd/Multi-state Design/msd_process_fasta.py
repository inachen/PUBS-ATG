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

# patches
PATCH1 = [1,19,20,57,60,62,63] 

PATCH2 = [42,44,45,46,48,49]

PATCH3 = [6,8,42,44,66,68,70,73]

PATCH4 = [44, 45, 57, 58, 59, 60, 62, 65] 

PATCH5 = [8, 42, 44, 68, 70, 72] #?

PATCHES = [PATCH1, PATCH2, PATCH3, PATCH4, PATCH5]

# UBQ
UBIQ_SEQ = "MQIFVKTLTGKTITLEVESSDTIDNVKSKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG!"
AA_LST= ["STOP", "W", "F", "Y", "L", "I", "M", "V", "C", "A", "G", "P", "S", "T", "N", "Q", "H", "R", "K", "D", "E"]

def get_fasta(filename):
    # converts file to fasta files
    infile = open(filename, 'r')
    outfilename = filename + ".fasta"
    outfile = open(outfilename, 'w')
    count = 0
    for line in infile:
        if ": AA:" in line and count < 125:
            count +=1
            outstring = str(">%s_%d\n") % (filename, count)
            fields = line.split()
            for element in fields:
                if "AA:" in element:
                    aa = element[-1]
                    outstring = outstring + aa
            outstring = outstring +"\n"
            outfile.write(outstring)
    outfile.close()

def get_prop_dic(fasta_fp, patch_num, run_num):

    prop_dic = {}

    patch_res = PATCHES[patch_num-1]

    # list of aa in each position of ubiquitin
    seq_lst = [[] for i in range(len(UBIQ_SEQ))]

    # list of sequences in string form
    seq_str_lst = []

    # storing strings for consensus
    seq_allel_lst = []

    # convert list of patch seq to ubiquitin seq list
    for record in SeqIO.parse(open(fasta_fp, "rU"), "fasta"):

        patch_seq = list(str(record.seq))
        ubi_seq = list(UBIQ_SEQ)
        allel_seq = ['-' for i in range(len(UBIQ_SEQ))]

        if len(patch_res) != len(patch_seq):
            raise ValueError("Patches don't match up")

        for pos, sub in zip(patch_res, patch_seq):

            ubi_seq[pos-1] = sub
            allel_seq[pos-1] = sub

        seq_str_lst.append(''.join(ubi_seq)[:-1])
        seq_allel_lst.append(''.join(allel_seq)[:-1])
        
        for i, aa in enumerate(ubi_seq):
            seq_lst[i].append(aa)


    # write WT sequences to list
    cw = csv.writer(open('out_txt/consensus_run_wt_'+str(run_num)+'_patch_'+str(patch_num)+'.csv', "wb"))

    for s in seq_str_lst:
        cw.writerow([s])

    # write sequences with only mutations
    cw = csv.writer(open('out_txt/consensus_run_allel_'+str(run_num)+'_patch_'+str(patch_num)+'.csv', "wb"))

    for s in seq_allel_lst:
        cw.writerow([s])

    # dictionary for consensus proportions
    prop_dic = {}

    # convert to dictionary 
    for pos, aa_lst in enumerate(seq_lst):

        total = len(aa_lst)
        for aa in AA_LST[1:]:
            num = aa_lst.count(aa)

            k = str(pos+1)
            prop_dic[k] = prop_dic.get(k, []) + ([(aa, num/float(total))])

    # write prop dic to pickle
    pickle.dump(prop_dic, open('outpickles/run_'+str(run_num)+'_patch_'+str(patch_num)+'_propdic.pkl', 'wb'))

def run():

    # file name format: run_#_patch_#_output
    filelst = glob("opt_output/run*_output.fasta")

    # goes through all output files
    for fp in filelst:

        info_lst = fp.split('_')
        run_num = info_lst[2]
        patch_num = int(info_lst[4])

        # # extract .fasta files from output
        # get_fasta(fp)

        # produce prop dictionary pickles (for bruk)
        # also saves sequence txt files for consensus logo
        get_prop_dic(fp, patch_num, run_num)

    # fasta_fp = "UQ_CON_interface_Patch1.fasta"

    # patch_num = 1

    # get_prop_dic(fasta_fp, patch_num)

run()


