import cPickle as pic
from enrichment_value import enrich_val
import numpy as np
num_dic = pic.load(open("aminotonumber.pkl", "rb"))

testset= pic.load(open("filtered_seq_perfect.pkl", "rb"))
wt_counts= pic.load(open("all_wt_barcode_counts.pkl", "rb"))

testset= enrich_val(testset)
########################## normalizing read enrichment by wt enrichment #################
ubiquitin = "MQIFVKTLTGKTITLEVESSDTIDNVKSKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
wt_aa_dict = dict(((i, ubiquitin[i-1]),[]) for i in range (1,77))
wt_aa_dict[(77, 'STOP')]= []
wt_aa_dict[(0, 'WT')]= []

for key in testset:
		for i in range (0,2):
			for j in range (0,3):
				testset[key][i][j] /= wt_counts[('enrich_wt')][i][j]
				
print testset[(4, 'F')]
				
for key in testset:
	if key in wt_aa_dict.keys():
		testset[key] = ([-1,-1,-1],[-1,-1,-1])
		
print testset[(4, 'F')]

############# dumping data into matrices ###################
i=int(raw_input("day: "))-1
j= int(raw_input("sample#: "))-1


fitnessarray = np.zeros((22,78))

print fitnessarray

for key in num_dic:
	print str(key) + ":" + str(num_dic[key])

for key in testset.iterkeys():
	if key[1] !='WT':
		aminonum = num_dic[key[1]]
	else:
		aminonum = 21
	fitnessarray[int(aminonum),int(key[0])] = testset[key][i][j]
print fitnessarray

filename = str(raw_input("file name: "))
# create a file called array file, write format
arrayfile= open(filename, "w")

#write to arrayfile from our array data with tab separators and numerical format
#fitnessarray.tofile(arrayfile, sep = "\t", format = "%F")

np.savetxt(arrayfile, fitnessarray, delimiter = "	")
pic.dump(fitnessarray, open("D2S3fitness_scores.py", "wb"))

