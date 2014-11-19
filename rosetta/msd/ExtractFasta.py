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

# patches

PATCH1 = [1,19,20,57,60,62,63] 

PATCH2 = [42,44,45,46,48,49]

PATCH3 = [6,8,42,44,66,68,70,73]

PATCH4 = [44, 45, 57, 58, 59, 60, 62, 65] 

PATCHES = [PATCH1, PATCH2, PATCH3, PATCH4]

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

def get_prop_dic(infasta):

    prop_dic = {}

    for record in SeqIO.parse(open(infasta, "rU"), "fasta"):
        print record



def run():

    filename = "UQ_CON_interface_Patch1.fasta"
    patch_num = 1

    get_prop_dic()


