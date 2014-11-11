import numpy as np
import cPickle as pic
data = pic.load(open("allele_dic_with_WT.pkl", "rb"))
alleles = []
alleles = data.keys()

def hamming_distance(s1, s2):
    """Calculate the Hammings Distance for equal length strings"""
    #Return the Hamming distance between equal-length sequences
    if len(s1) != len(s2):
        # raise ValueError("Undefined for sequences of unequal length")
        return -1

    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))
    
#comparing two lists of hamming distances
'''def mult_comp_list(seq1, seq2):
	#hamming distance of 1 or 2 that encode for the same or different mutations
	num_ham_same = 0
	num_ham_diff = 0
	for i in range(0, len(seq1)):
		for j in range(0, len(seq2)):
			n = 0
			n = hamming_distance(seq1[i], seq2[j])
			if n > 0 and n < 3:
				if se
			#num_list.append(n)
	return num_list'''
	
print data[2].key