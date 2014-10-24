import cPickle as pic
from enrichment_value import enrich_val
import numpy as np
num_dic = pic.load(open("aminotonumber.pkl", "rb"))

testset= pic.load(open("filtered_seq_perfect.pkl", "rb"))
wt_counts= pic.load(open("wt_counts.pkl", "rb"))
testset= enrich_val(testset)

for key in testset:
		for i in range (0,2):
			for j in range (0,3):
				testset[key][i][j] /= wt_counts[('enrich_wt')][i][j]
		
pic.dump(testset, open("fitness_scores.pkl", "wb"))	

print testset


i=int(raw_input("day: "))-1
j= int(raw_input("sample#: "))-1


fitnessarray = np.zeros((21,78))
for key in testset.iterkeys():
	if key[1] !='WT':
		aminonum = num_dic[key[1]]
	else:
		aminonum = 21
	fitnessarray[int(aminonum-1),int(key[0])] += testset[key][i][j]
print fitnessarray

filename = str(raw_input("file name: "))
# create a file called array file, write format
arrayfile= open(filename, "w")

#write to arrayfile from our array data with tab separators and numerical format
#fitnessarray.tofile(arrayfile, sep = "\t", format = "%F")

np.savetxt(filename, fitnessarray, delimiter = "	")

