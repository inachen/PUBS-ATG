import numpy as np
import matplotlib
import cPickle as pic
seq_data = pic.load(open("2500_filtered_seq.pkl", "rb"))

'''#making a dictionary called wildtypedict to identify wt
ubiquitin = "MQIFVKTLTGKTITLEVESSDTIDNVKSKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
#testset = pic.load(open("2500_filtered_seq.pkl", "rb"))
wildtypedict = dict(((i, ubiquitin[i-1]),([],[])) for i in range (1,76))
wildtypedict[(76, 'STOP')]=([],[])
for key in wildtypedict.keys():
	if key in seq_data.keys():
		wildtypedict[key] = seq_data[key]'''



#takes dictionary and calculates total counts at each position
def total_value(seq_data):
	sum_counts = {('total') : ([0,0,0],[0,0,0])}
	for keys in seq_data.iterkeys():
		for i in range(0,1):
			for j in range(0,2):
				sum_counts['total'][i][j] += seq_data[keys][i][j]
	return sum_counts
				
#takes in a dictionary, form: ((key):([N,N,N],[N,N,N])) and returns Enrichment values at each time point in the same format
def enrich_val(seq_data):
	#sum_counts is the total counts in a dictionary
	sum_counts = total_value(seq_data)
	enrich_val = {}
	#iterating over all the values in all the keys and doing the enrichment value calculations (returns as a float)
	for keys in seq_data.iterkeys():
		#declaring dictionary to use for values
		enrich_val[keys] = ([0,0,0],[0,0,0])
		for i in range(0,1):
			for j in range(0,2):
				F = 0.0
				F_zero = 0.0
				E = 0.0
				F_zero = float(seq_data[keys][i][0])/float(sum_counts['total'][i][0])
				F = float(seq_data[keys][i][j])/float(sum_counts['total'][i][j])
				#need to remove zero value so error for dividing by zero does not occur
				if F_zero != 0:
					E = float(F/F_zero)
					enrich_val[keys][i][j] = E
	return enrich_val
		

#wildtype_E = enrich_val(wildtypedict)
#print enrich_val(seq_data)
#print wildtype_E
#def enrich_val(seq_data):
	