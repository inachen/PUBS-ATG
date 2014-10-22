import sys
from Bio import SeqIO
import cPickle as pic
from barcode_differences import single_comp

ATG_indices = ['TTGACT','GGAACT','TGACAT','GGACGG','CTCTAC','GCGGAC']
filename = raw_input("sequencing fastq file name: ")
filename = "testseq.fastq"

count = 0
indexcount = 0
theoutput= dict((element,[]) for element in ATG_indices )
indexcounter = dict((element,[]) for element in ATG_indices)
#print theoutput
#print indexcounter

for value in ATG_indices:
	for record in SeqIO.parse(filename, "fastq"):
		indexread = record.description.split(":")[-1]
		indexcounter[value].append((indexread,0))
		if single_comp(value, indexread) <1:
			theoutput[value].append((record.seq,record.letter_annotations["phred_quality"]))
			#indexcounter[value]]+=1
			
#print indexcounter
pic.dump(theoutput, open("theoutput.pkl", "wb"))