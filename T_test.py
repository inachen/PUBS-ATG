import cPickle as pic
import numpy as np

fit_scores = pic.load(open("fitness_scores.pkl", "rb"))

mean_1 = [0] * 80
mean_2 = []
num_counts_1 = []
num_counts_2 = []
stdev_1 = []
stdev_2 = []
s_pooled = []

aa_num = 0
for keys in fit_scores.iterkeys():
	for other_keys in fit_scores.iterkeys():
		if other_keys[0] == keys[0]:
			mean_1[keys[0]-1] += fit_scores[keys][0][2]
			mean_1 = np.mean(mean_1)

def average_cal(num_list):
	mean = 0
	np.mean(num_list)
	return mean

def stdev_calc(num_list):
	standard_dev = 0
	np.std
	
#def S_pooled(full_dict):
	
print mean_1
#print fit_scores