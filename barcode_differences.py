barcodes = ['CGTGAT', 'ACATCG', 'GCCTAA','TGGTCA','CACTGT','ATTGGC','GATCTG','TCAAGT','CTGATC','AAGCTA','GTAGCC','TACAAG','TTGACT','GGAACT','TGACAT','GGACGG','CTCTAC','GCGGAC','TTTCAC','GGCCAC','CGAAAC','CGTACG','CCACTC','GCTACC','ATCAGT','GCTCAT','AGGAAT','CTTTTG','TAGTTG','CCGGTG']
import numpy as np
import cPickle as pic
data = pic.load(open("allele_dic_with_WT.pkl", "rb"))
alleles = []
alleles = data.keys()

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

#looking through two lists of strings, and two that have > 2 differences will be the output
def check(seq1, seq2):
	mat = mult_comp(seq1, seq2)
	f = []
	n = 0
	for i in range (0, len(seq1)):
		for j in range (0, len(seq2)):
			if mat[i,j] > 2:
				f.append(str(i + 1) + ": " + seq1[i] + ", " +  str(j + 1) +": " + seq2[j])
	return f			
	
#print single_comp("cats", "dogs")

#print mult_comp(alleles, alleles)
#print check(barcodes, barcodes)
