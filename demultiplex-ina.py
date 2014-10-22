import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import collections
import csv
import numpy as np
import math
import matplotlib.cm as cm
from mpltools import style
from mpltools import color
from Bio import SeqIO
import cPickle as pickle
# from barcode_differences.py import single_comp

# Constants
ATG_INDICES = ['TTGACT','GGAACT','TGACAT','GGACGG','CTCTAC','GCGGAC']

# BARCODES = 

def demultiplex ():
	infile = "testseq.fastq"
	# newpickle = pic.load(open("newpickle.pkl", "rb"))

	# create a dictionary to store the barcodes
	seqdata = dict((i,{}) for i in ATG_INDICES)

	for record in SeqIO.parse(open(infile, "rU"), "fastq"):
		
		# pull out index from description
		index = record.description.split(":")[-1]
		sequence = str(record.seq)

		# update counts
		if index in seqdata:
			seqdata[index][sequence] = seqdata[index].get(sequence, 0) + 1


	outfile = open('rawseq.txt', 'wb')
	pickle.dump(seqdata, outfile)
	# return seqdata

def barcode_filter():

	infile = open('rawseq.txt', 'rb')
	seqdata = pickle.load(infile)

	print seqdata
			
def run():
	demultiplex()
	barcode_filter()

run()

