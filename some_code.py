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
from Bio import SeqIO
import cPickle as pickle

############ TO DO ################
# 1. add near match to indices
# 2. bar_to_aa()
# 3. make_aa_dic()
###################################

############ Functions ############
# 
###################################

# bad illumina quality threshold
BAD_QUAL_THRESHOLD = 20

# maximum number of bad quality reads in barcode
BAD_QUAL_NUM = 3

# length of barcode sequence
BARCODE_LEN = 18

# hammings distance error tolerances 
BARCODE_ERR = 1
INDEX_ERR = 1

# ambiguity threshold (other sequences farther than this away)
BARCODE_SEP = 4

# namedtuple to represent timepoint objects
# day = [1,2], time = [1,2,3], index = sequence
Timepoint = namedtuple("Timepoint", ["day", "time", "seq"])

# Constants
ATG_INDICES = [Timepoint(day=1, time=1, seq='TTGACT'),
            Timepoint(day=1, time=2, seq='GGAACT'),
            Timepoint(day=1, time=3, seq='TGACAT'),
            Timepoint(day=2, time=1, seq='GGACGG'),
            Timepoint(day=2, time=2, seq='CTCTAC'),
            Timepoint(day=2, time=3, seq='GCGGAC')]

# source: wikipedia
def hamming_distance(s1, s2):
    #Return the Hamming distance between equal-length sequences
    if len(s1) != len(s2):
        # raise ValueError("Undefined for sequences of unequal length")
        return -1

    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def bar_to_aa(bar):
    return 

def make_aa_dic():
    return

def demultiplex ():
    infile = "testseq.fastq"

    # create a dictionary to store the barcodes
    seqdata = dict((i,[]) for i in ATG_INDICES)
    index_lst = [i.seq for i in ATG_INDICES]

    for record in SeqIO.parse(open(infile, "rU"), "fastq"):
        
        # pull out index from description
        index_seq = record.description.split(":")[-1]

        # get reverse complement of the barcode!
        sequence = str(record.seq.reverse_complement())[7:]
        qual = record.letter_annotations["phred_quality"]

        # update counts
        if index_seq in index_lst:
            # seqdata[index][sequence] = seqdata[index].get(sequence, 0) + 1
            seqdata[ATG_INDICES[index_lst.index(index_seq)]].append((sequence, qual))

    # dump dictionary into a pickle
    outfile = open('rawseq.pkl', 'wb')
    pickle.dump(seqdata, outfile)

def barcode_filter():

    # open seq data
    infile = open('rawseq.pkl', 'rb')
    seqdata = pickle.load(infile)

    # get barcode library
    bar_dic = pickle.load(open("allele_dic_with_WT.pkl", "rb"))
    bar_lst = bar_dic.keys()

    # trajectory dictionary
    bar_timecourse = dict((i, ([0,0,0], [0,0,0])) for i in bar_lst)

    # counts different types of matches
    count_match = 0
    count_near = 0
    count_ambig = 0
            
    # update timecourse dictionary
    for timept, seqlst in seqdata.iteritems():

        for s in seqlst:
            sequence = s[0]
            qual = s[1]

            # check for bad quality
            if sum(i < BAD_QUAL_THRESHOLD for i in qual) > BAD_QUAL_NUM :
                continue

            # check for length
            if len(sequence) != BARCODE_LEN :
                continue

            # store match sequence
            match = ""
            near_matches = []

            # check for match
            for bar in bar_lst:

                d = hamming_distance(bar, sequence)

                if d == 0:
                    match = bar
                    break

                elif d <= BARCODE_ERR + BARCODE_SEP:
                    near_matches.append((bar, d))

            # update match

            # exact match
            if match != "":

                bar_timecourse[match][timept.day - 1][timept.time - 1] += 1

                count_match +=1

            # one close match within error range
            elif (len(near_matches) == 1) and (near_matches[0][1] <= BARCODE_ERR):
                match_seq = near_matches[0][0]
                bar_timecourse[match_seq][timept.day - 1][timept.time - 1] += 1

                count_near += 1

            # more than one matches
            elif len(near_matches) > 1:

                # sort by distance
                sorted_matches = sorted(near_matches, key=lambda x: x[1])
                first_dist = sorted_matches[0][1]
                second_dist = sorted_matches[1][1]

                # first within range, second farther away than ambiguity threshold
                if (first_dist <= BARCODE_ERR) and (second_dist - first_dist > BARCODE_SEP) :
                    match_seq = sorted_matches[0][0]
                    bar_timecourse[match_seq][timept.day - 1][timept.time - 1] += 1

                    count_near += 1

                else :
                    count_ambig += 1

    # print bar_timecourse

    print "Identical matches: ", count_match
    print "Near matches: ", count_near
    print "Ambiguous cases: ", count_ambig
            
def run():
    # demultiplex()
    barcode_filter()

run()

