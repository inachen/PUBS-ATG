
######## descriptions BELOW code ########

"""this program takes in a fastq file, demultiplexes according to provided indices, 
demultiplexes with a distance treshold of 1, outputs a pickled demultiplexed dictionary 
("demultiplexed_reads.pkl")and prints a summary of the demultiplexing"""


import sys
from Bio import SeqIO
import cPickle as pic
from barcode_differences import single_comp
from barcode_differences import mult_comp
import numpy as np

ATG_indices = ['TTGACT','GGAACT','TGACAT','GGACGG','CTCTAC','GCGGAC']

filename = raw_input("sequencing fastq file name: ")
#asks for the name of the fastq file you want to demultiplex

theoutput= dict((element,[]) for element in ATG_indices )
#theoutput = the dictionary of demultiplexed reads with 2plets of (sequence, processed quality string)

indexcounter = dict((element,dict()) for element in ATG_indices)
#for counting the number of perfect and non-perfect matches that binned into each index

total_reads= dict((element,0) for element in ATG_indices)
#for counting the total reads for each time-point index

for value in ATG_indices:
#iterates over the ATG_indices

	for record in SeqIO.parse(filename, "fastq"):
	#iterates over lines of sequence in input fastq file
	
		indexread = record.description.split(":")[-1]
		# isolates the index from the description line of sequence
		
		indexcounter[value][indexread]=indexcounter[value].get(indexread,1)+1
		# creates an entry in indexcounter if index is new, gives it a 1 count, otherwise finds the index and increments the count
		
		if single_comp(value, indexread) <1:
		# checks if the index from seq is 1 distance or less from the desired index
		
			theoutput[value].append((record.seq,record.letter_annotations["phred_quality"]))
			#appends the sequence and quality of an approved read into the desired index bin
			
			total_reads[value]+= 1
			# adds to the total read count for the relevant index
			
pic.dump(theoutput, open("demultiplexed_reads.pkl", "wb"))
#saves the output dictionary as a pickle

print total_reads