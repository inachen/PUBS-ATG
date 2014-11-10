import cPickle as pic
import numpy as np

wt_aa_dict = pic.load(open("wt_aa_dict.pkl", "rb"))
wt_number_dict = pic.load(open("wt_number_dict.pkl", "rb"))
testset = pic.load(open("test_set.pkl", "rb"))

print "printing testset (57, 'D') before norm"
print testset[57, 'D']
print "printing testset (57, 'S') before norm"
print testset[57, 'S']
print "normalization factor"
print wt_number_dict[57]

for key in testset:
	if key[0] != 0:
		aa_pos = key[0]	
		for i in range (0,2):
			for j in range (0,3):
				#print aa_pos
				testset[key][i][j] /= wt_number_dict[aa_pos][i][j]
		
pic.dump(wt_aa_dict, open("fitness_scores", "wb"))	
# adds wt counts for aa position to mutant counts for aa position in testset

print "printing testset (57, 'D') after norm"
print testset[57, 'D']
print "printing testset (57, 'S') after norm"
print testset[57, 'S']
