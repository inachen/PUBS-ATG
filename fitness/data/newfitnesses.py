import numpy as np
import cPickle as pic
num_dic = pic.load(open("aminotonumber.pkl", "rb"))
aa_dic = pic.load(open("translate.pkl", "rb"))

#################### change name of file to input ##########################
testset = pic.load(open("Codon_FitScore_MG132_day1.pkl", "rb"))

fitnessarray = np.zeros((21,77))


for key in testset.iterkeys():
    rnacodon = key[1].replace('T','U')
    if rnacodon != 'WU':
		aa_iden = aa_dic[rnacodon]
		aminonum = num_dic[aa_iden]
		fitnessarray[int(aminonum-1),int(key[0]-1)] = testset[key][0]

################### Change the name of the output file ######################
filename = "MG132_D1_new_fitness_scores.txt"

arrayfile= open(filename, "w")

################### Change the name of the output pickle ###################
np.savetxt(arrayfile, fitnessarray, delimiter = ",")
pic.dump(fitnessarray, open("MG132_D1_new_fitness_scores.pkl", "wb"))