import numpy as np
import cPickle as pic
fitness_dict = pic.load(open("input/D2S3fitness_scores.pkl", "rb"))
ddG_dict = pic.load(open("out_pickles/ddg_dic.pkl", "rb"))


hydrophobic = ['M', 'C', 'I', 'L', 'Y', 'F', 'W']
polar = ['Q', 'P', 'N', 'A', 'T', 'S', 'V', 'G']
charged = ['E', 'D', 'K', 'R', 'H']
compiled = {}

#print fitness_dict
count = 0
for key in ddG_dict:
	try:
		fitness_score = fitness_dict[key]
		ddG_score = ddG_dict[key]
		compiled[key] = (fitness_score, ddG_score)
	except KeyError as error:
		count +=1
		print error
		print ddG_dict[key]
		print key in fitness_dict.keys()
		pass
	
	
print count
	
	