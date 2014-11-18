import numpy as np
import cPickle as pic
num_dic = pic.load(open("aminotonumber.pkl", "rb"))
aa_dic = pic.load(open("translate.pkl", "rb"))
testset = pic.load(open("Codon_FitScore_HU_day1.pkl", "rb"))

fitnessarray = np.zeros((22,78))


print fitnessarray

for key in testset.iterkeys():
	codon = key[1]
	rnacodon = codon.replace('T','U')
	if rnacodon != 'WU':
		aa_iden = aa_dic[rnacodon]
	else:
		aminonum = 21
	if key[1] !='WT':
		aminonum = num_dic[aa_iden]
		fitnessarray[int(aminonum),int(key[0])] = testset[key][0]
print fitnessarray[1]

filename = "D1_new_fitness_scores.txt"
# create a file called array file, write format
arrayfile= open(filename, "w")

#write to arrayfile from our array data with tab separators and numerical format
#fitnessarray.tofile(arrayfile, sep = "\t", format = "%F")

np.savetxt(arrayfile, fitnessarray, delimiter = "	")
pic.dump(fitnessarray, open("D1_new_fitness_scores.py", "wb"))