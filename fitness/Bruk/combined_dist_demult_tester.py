import numpy as np
import matplotlib.pyplot as plt
import cPickle as pic
data = pic.load(open("allele_dic_with_WT.pkl", "rb"))
alleles = []
alleles = data.keys()

import sys
from Bio import SeqIO

#comparing one string to another and returning the number of differences
def single_comp(cats, dogs):
	n = 0
	if cats == dogs:
		return n
	else:
		for i in range(0, len(cats)):
			if cats[i] != dogs[i]:
				n += 1
		return n

#comparing two lists of strings and putting the number of differences into a comparative matrix
def mult_comp(seq1, seq2):
	mat = (len(seq1), len(seq2))
	mat = np.zeros(mat)
	for i in range(0, len(seq1)):
		for j in range(0, len(seq2)):
			n = 0
			n = single_comp(seq1[i], seq2[j])
			mat[i,j] = n
	return mat
	
ATG_indices = ['TTGACT','GGAACT','TGACAT','GGACGG','CTCTAC','GCGGAC']

thefilename = ["100k_testdata.fastq"]
#this line is for testing on small datasets, ignore!

#thefilename = ["/data/ClassData-2014/data-oct16/lane1_Undetermined_L001_R1_001.fastq", "/data/ClassData-2014/data-oct16/lane2_Undetermined_L002_R1_001.fastq", "/data/ClassData-2014/data-oct17/lane1_Undetermined_L001_R1_001.fastq"]
#files to demultiplex

#filename = raw_input("sequencing fastq file name: ")
#optional manual input for file names, ignore!

theoutput= dict((element,[]) for element in ATG_indices )
#theoutput = the dictionary of demultiplexed reads with 2plets of (sequence, processed quality string)

indexcounter = dict((element,dict()) for element in ATG_indices)
#for counting the number of perfect and non-perfect matches that binned into each index

binindexcounter = dict((element,dict()) for element in ATG_indices)

total_reads= dict((element,0) for element in ATG_indices)
#for counting the total reads for each time-point index

for filename in thefilename:
	print "processing" + filename
	for value in ATG_indices:
	#iterates over the ATG_indices

		for record in SeqIO.parse(filename, "fastq"):
		#iterates over lines of sequence in input fastq file
	
			indexread = record.description.split(":")[-1]
			# isolates the index from the description line of sequence
		
			indexcounter[value][indexread]=indexcounter[value].get(indexread,1)+1
			# creates an entry in indexcounter if index is new, gives it a 1 count, otherwise finds the index and increments the count
		
			if single_comp(value, indexread) <2:
			# checks if the index from seq is 1 distance or less from the desired index
		
				theoutput[value].append((record.seq,record.letter_annotations["phred_quality"]))
				#appends the sequence and quality of an approved read into the desired index bin
			
				total_reads[value]+= 1
				# adds to the total read count for the relevant index
				
				binindexcounter[value][indexread]=binindexcounter[value].get(indexread, 1)+1
				# creates an entry in binindexcounter if binned index is new, gives it a 1 count, otherwise find and increments index
				
			
pic.dump(theoutput, open("demultiplexed_reads.pkl", "wb"))
#saves the output dictionary as a pickle

print total_reads
print binindexcounter