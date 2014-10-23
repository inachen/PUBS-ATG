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

############ TO DO ################
# 1. add near match to indices
###################################

############ Functions ############
# 
###################################

# bio constants
NUM_AA = 21
NUM_POS = 77

# bad illumina quality threshold
BAD_QUAL_THRESHOLD = 20

# maximum number of bad quality reads in barcode
BAD_QUAL_NUM = 3

# length of barcode sequence
INDEX_LEN = 6
BARCODE_LEN = 18

# hammings distance error tolerances 
BARCODE_ERR = 1
INDEX_ERR = 1

# ambiguity threshold (other sequences farther than this away)
BARCODE_SEP = 3
INDEX_SEP = 3

# directory names
OUTPUT_DIR = "output/"
PICKLE_DIR = "pickles/"
TESTSET_DIR = "test_sets/"

# reference pickle files
ALLELE_PIC = "allele_dic_with_WT.pkl"
TRANSLATE_PIC = "translate.pkl"

# output pickle files
BAR_AA_PIC = "bar_to_aa.pkl"
DEMULT_PIC = "rawseq.pkl"
FILT_PIC = "filtered_seq.pkl"

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

FASTA_INPUT = ["/data/ClassData-2014/data-oct16/lane1_Undetermined_L001_R1_001.fastq", 
    "/data/ClassData-2014/data-oct16/lane2_Undetermined_L002_R1_001.fastq", 
    "/data/ClassData-2014/data-oct17/lane1_Undetermined_L001_R1_001.fastq"]

TEST_INPUT = [TESTSET_DIR + "2500_testdata.fastq"]

# source: wikipedia
def hamming_distance(s1, s2):
    #Return the Hamming distance between equal-length sequences
    if len(s1) != len(s2):
        # raise ValueError("Undefined for sequences of unequal length")
        return -1

    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def make_aa_dic():

    # reference dictionaries
    bar_codon_dic = pickle.load(open(PICKLE_DIR + ALLELE_PIC, "rb"))
    codon_aa_dic = pickle.load(open(PICKLE_DIR + TRANSLATE_PIC, "rb"))

    # make dictionary to store barcode to amino acid table
    bar_aa_dic = copy.deepcopy(bar_codon_dic)

    for bar, pair in bar_aa_dic.iteritems():
        if pair[1] != 'WT':
            bar_aa_dic[bar] = (pair[0], codon_aa_dic[pair[1].replace('T', 'U')])

    pickle.dump(bar_aa_dic, open(OUTPUT_DIR + BAR_AA_PIC, 'wb'))

    print "======= Made Barcode to AA Dictionary ======="

def demultiplex (fileset):

    print "====== Started Demultiplexing ======"

    # create a dictionary to store the barcodes
    seqdata = dict((i,[]) for i in ATG_INDICES)
    index_lst = [i.seq for i in ATG_INDICES]

    count_match = 0
    count_near = 0
    count_ambig = 0

    count_total = 0
    count_bad_qual = 0
    count_bad_len = 0

    for infile in fileset:
        print "Demultiplexing ", infile, "..."
        for record in SeqIO.parse(open(infile, "rU"), "fastq"):

            count_total += 1
            
            # pull out index from description
            index_seq = record.description.split(":")[-1]

            # get reverse complement of the barcode!
            sequence = str(record.seq.reverse_complement())[7:]
            qual = record.letter_annotations["phred_quality"]

            # check for bad quality
            if sum(i < BAD_QUAL_THRESHOLD for i in qual) > BAD_QUAL_NUM :
                count_bad_qual += 1
                continue

            # check for length
            if len(index_seq) != INDEX_LEN :
                count_bad_len += 1
                continue

            # store match sequence
            match = ""
            near_matches = []

            # check for match
            for ind in index_lst:

                d = hamming_distance(ind, index_seq)

                if d == 0:
                    match = ind
                    break

                elif d <= INDEX_ERR + INDEX_SEP:
                    near_matches.append((ind, d))

            # update match

            # exact match
            if match != "":

                seqdata[ATG_INDICES[index_lst.index(match)]].append((sequence, qual))

                count_match +=1

            # one close match within error range
            elif (len(near_matches) == 1) and (near_matches[0][1] <= INDEX_ERR):
                match_ind = near_matches[0][0]
                seqdata[ATG_INDICES[index_lst.index(match_ind)]].append((sequence, qual))

                count_near += 1

            # more than one matches
            elif len(near_matches) > 1:

                # sort by distance
                sorted_matches = sorted(near_matches, key=lambda x: x[1])
                first_dist = sorted_matches[0][1]
                second_dist = sorted_matches[1][1]

                # first within range, second farther away than ambiguity threshold
                if (first_dist <= INDEX_ERR) and (second_dist - first_dist > INDEX_SEP) :
                    match_ind = sorted_matches[0][0]
                    seqdata[ATG_INDICES[index_lst.index(match_ind)]].append((sequence, qual))

                    count_near += 1

                else :
                    count_ambig += 1

            # # update counts
            # if index_seq in index_lst:
            #     # seqdata[index][sequence] = seqdata[index].get(sequence, 0) + 1
            #     seqdata[ATG_INDICES[index_lst.index(index_seq)]].append((sequence, qual))

    # dump dictionary into a pickle
    pickle.dump(seqdata, open(OUTPUT_DIR + DEMULT_PIC, 'wb'))

    print"======= Finished Demultiplexing ======="

    print "====== Index Filter Counts ======"
    print "Total Count :", count_total
    print "Bad Quality Count: ", count_bad_qual
    print "Bad Length Count: ", count_bad_len
    print "Identical Matches :", count_match
    print "Near Matches :", count_near
    print "Ambiguous Cases: ", count_ambig

def check_demult():

    print "====== Demultiplexing Counts ======"

    infile = open(OUTPUT_DIR + DEMULT_PIC, 'rb')
    seqdata = pickle.load(infile)
    infile.close()

    # checks how many counts there are for each index
    for index, seqlst in seqdata.iteritems():
        print index.seq, ': ', len(seqlst)

def barcode_filter():

    print "======= Started Barcode Filtering ======="

    # open seq data
    infile = open(OUTPUT_DIR + DEMULT_PIC, 'rb')
    seqdata = pickle.load(infile)
    infile.close()

    # get barcode lst
    bar_dic = pickle.load(open(PICKLE_DIR + ALLELE_PIC, "rb"))
    bar_lst = bar_dic.keys()

    # get list of unique mutations
    bar_aa_dic = pickle.load(open(OUTPUT_DIR + BAR_AA_PIC, 'rb'))
    mut_lst = list(set(bar_aa_dic.values()))

    # trajectory dictionary 
    # ([timept 1, timept 2, timept 3], [timept 1, timept 2, timept 3])
    # first list is day 1, second list is day two
    mut_timecourse = dict((i, ([0,0,0], [0,0,0])) for i in mut_lst)

    # counts different types of matches
    count_match = 0
    count_near = 0
    count_ambig = 0

    count_total = 0
    count_bad_qual = 0
    count_bad_len = 0
            
    # update timecourse dictionary
    for timept, seqlst in seqdata.iteritems():

        print "Filtering", timept, "..."

        for s in seqlst:
            sequence = s[0]
            qual = s[1]

            count_total += 1

            # check for bad quality
            if sum(i < BAD_QUAL_THRESHOLD for i in qual) > BAD_QUAL_NUM :
                count_bad_qual += 1
                continue

            # check for length
            if len(sequence) != BARCODE_LEN :
                count_bad_len += 1
                continue

            # store match sequence
            match = ""
            near_matches = []

            # check for match
            for bar in bar_lst:

                d = hamming_distance(bar, sequence)

                if d == 0:
                    match = bar_aa_dic[bar]
                    break

                elif d <= BARCODE_ERR + BARCODE_SEP:
                    near_matches.append((bar_aa_dic[bar], d))

            # update match

            # exact match
            if match != "":

                mut_timecourse[match][timept.day - 1][timept.time - 1] += 1

                count_match +=1

            # one close match within error range
            elif (len(near_matches) == 1) and (near_matches[0][1] <= BARCODE_ERR):
                match_seq = near_matches[0][0]
                mut_timecourse[match_seq][timept.day - 1][timept.time - 1] += 1

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
                    mut_timecourse[match_seq][timept.day - 1][timept.time - 1] += 1

                    count_near += 1

                else :
                    count_ambig += 1

    print mut_timecourse

    pickle.dump(mut_timecourse, open(OUTPUT_DIR + FILT_PIC, 'wb'))

    print "====== Completed Barcode Filtering ======"
    print "====== Barcode Filter Counts ======"
    print "Total Count :", count_total
    print "Bad Quality Count: ", count_bad_qual
    print "Bad Length Count: ", count_bad_len
    print "Identical Matches :", count_match
    print "Near Matches :", count_near
    print "Ambiguous Cases: ", count_ambig

            
def run():

    ### RUN ONCE - makes rawseq.pkl ###
    if not (os.path.exists(OUTPUT_DIR + DEMULT_PIC)):

        # FASTA_INPUT for actual data, TEST_INPUT for small test set
        demultiplex(fileset = TEST_INPUT)
        check_demult()

    ### RUN ONCE - makes bar_to_aa.pkl ###
    if not (os.path.exists(OUTPUT_DIR + BAR_AA_PIC)):
        make_aa_dic()

    ### RUN ONCE - make the amino acid table of filtered counts ###
    if not (os.path.exists(OUTPUT_DIR + FILT_PIC)):
        barcode_filter()


run()

