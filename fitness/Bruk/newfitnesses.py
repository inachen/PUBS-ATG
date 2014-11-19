import numpy as np
import cPickle as pic
num_dic = pic.load(open("aminotonumber.pkl", "rb"))
aa_dic = pic.load(open("translate.pkl", "rb"))

#################### change name of file to input ##########################
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

####
filename = "HU_D1_new_fitness_scores.txt"

arrayfile= open(filename, "w")


np.savetxt(arrayfile, fitnessarray, delimiter = "	")
pic.dump(fitnessarray, open("HU_D1_new_fitness_scores.py", "wb"))