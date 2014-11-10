import cPickle as pic
from enrichment_value import enrich_val
import numpy as np
num_dic = pic.load(open("aminotonumber.pkl", "rb"))

testset= pic.load(open("filtered_seq_perfect.pkl", "rb"))
num_dic = pic.load(open("aminotonumber.pkl", "rb"))

total_counts = {'total':([0.0,0.0,0.0],[0.0,0.0,0.0])}

#add up all counts
for key in testset:
	for i in range (0,2):
		for j in range (0,3):
			total_counts['total'][i][j] += testset[key][i][j]

#print total_counts

####################### for fractional count matrix ##############################

fractional_count = {}
print fractional_count


for key in testset:
	fractional_count[key]=([0.0,0.0,0.0],[0.0,0.0,0.0])
	for i in range (0,2):
		for j in range (0,3):
			fractional_count[key][i][j] = testset[key][i][j]/total_counts['total'][i][j]
#print fractional_count

i=int(raw_input("day: "))-1
j= int(raw_input("sample#: "))-1


frac_array = np.zeros((22,78))

print frac_array

for key in fractional_count.iterkeys():
	if key[1] !='WT':
		aminonum = num_dic[key[1]]
	else:
		aminonum = 21
	frac_array[int(aminonum),int(key[0])] = fractional_count[key][i][j]
print frac_array

print testset[4, 'F']
print total_counts
print fractional_count[4, 'F']
print frac_array[2,4]

filename = str(raw_input("file name: "))
# create a file called array file, write format
arrayfile= open(filename, "w")

#write to arrayfile from our array data with tab separators and numerical format
#fitnessarray.tofile(arrayfile, sep = "\t", format = "%F")

np.savetxt(arrayfile, frac_array, delimiter = "	")
#####################################################################################