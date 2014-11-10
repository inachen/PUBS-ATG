import cPickle as pic
import numpy as np
bar_nuc = pic.load(open("allele_dic.pkl", "rb"))
codon_dic = pic.load(open("translate.pkl", "rb"))
num_dic = pic.load(open("aminotonumber.pkl", "rb"))

count = 0
freqarray = np.zeros((20,77))
for key in bar_nuc.iterkeys():
	# call up value for key, split by _, assign split values back to keys
	bar_nuc[key] = bar_nuc[key][0].split("_")
	# pull out second half of split, name it dnacodon
	dnacodon= bar_nuc[key][1]
	# replaces T's with U's in dnacodon, names it rnacodon
	rnacodon = dnacodon.replace('T','U')
	# calls up the values associated with 'rnacodon' key, names it amino
	amino = codon_dic[rnacodon]
	aminonum = num_dic[amino]
	# replaces the 2nd value of each key with amino
	bar_nuc[key][1] = aminonum
	#incrementing each incidence
	freqarray[int(bar_nuc[key][1]-1),int(bar_nuc[key][0])-1]+=1
	
	
	
	
print array
# create a file called array file, write format
arrayfile= open("arrayfile", "w")

#write to arrayfile from our array data with tab separators and numerical format
freqarray.tofile(arrayfile, sep = "\t", format = "%d")


