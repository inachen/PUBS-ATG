import numpy as np
import cPickle as pic
fitness_dict = pic.load(open("D2S3fitness_scores.pkl", "rb"))
ddG_dict = pic.load(open("ddg_dic.pkl", "rb"))


hydrophobic = ['M', 'C', 'I', 'L', 'Y', 'F', 'W']
polar = ['Q', 'P', 'N', 'A', 'T', 'S', 'V', 'G']
charged = ['E', 'D', 'K', 'R', 'H']
compiled = {}

#print fitness_dict

for key in fitness_dict:
	if key[1] != 'STOP':
		fitness_score = fitness_dict[key]
		ddG_score = ddG_dict[key]
		compiled[key] = (fitness_score, ddG_score)
		print compiled[key]
	
	
print compiled
	
	