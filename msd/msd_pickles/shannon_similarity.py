import cPickle as pic
import numpy as np
from glob import glob
import math
from collections import Counter

patch_1 = ['1','19','20','57','60','62','63']
patch_2 = []
patch_3 = []
patch_4 = []

patch = patch_1

def entropy(s):
	p, lns = Counter(s), float(len(s))
	return -sum( count/lns * math.log(count/lns, 2) for count in p.values())

filelist = glob("run*")

distance_matrix = np.zeros((len(filelist),len(filelist)))
for i in range(len(filelist)):
	for j in range(len(filelist)):
		rosetta_logo = pic.load(open(filelist[i],"rb"))
		fitness_logo = pic.load(open(filelist[j],"rb"))
		rosetta_list = []
		fitness_list = []
		difference = []
		distance_array = []
		average_array = []
		final_score = 0
		patch1 = filelist[i][12]
		patch2
		if i != j:
			print i,j
			for key in rosetta_logo:
				if key in patch:
					for r, f in zip(rosetta_logo[key], fitness_logo[key]):
						rosetta_list.append(r)
						fitness_list.append(f)
					#print rosetta_list
					#print fitness_list

					sorted(rosetta_list)
					sorted(fitness_list)
					#print rosetta_list
					#print fitness_list

			
					for k in range(len(rosetta_list)):
						difference.append(rosetta_list[k][1]-fitness_list[k][1])
						distance_array.append(entropy(difference))
					
					average_dist = sum(distance_array)/(len(distance_array)+1)
					average_array.append(average_dist)
			print average_array
			final_score = sum(average_array)/len(average_array)
			print "final score", final_score
			distance_matrix[i,j] = final_score

print distance_matrix





