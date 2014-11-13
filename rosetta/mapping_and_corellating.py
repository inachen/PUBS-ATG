fake_list = ([1,2,3,1,2,3,1,2,4],[1,2,4,5,7,9,2,3,4])
import csv
import numpy as np
import matplotlib.pyplot as plt
fitness_lists = open("out_pickles/all_DDG_Fitness.csv", "r")

def correlation(data):
	final = np.corrcoef(data[0], data[1])[0][1]
	return final

def dot_plot(data):
	plt.plot(data[0], data[1], "ro")
	plt.ylabel('some numbers')
	plt.show()


#print correlation(fake_list)
#dot_plot(fake_list)
print fitness_lists

