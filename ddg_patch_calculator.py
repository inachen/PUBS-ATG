import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import collections
from collections import namedtuple
import csv
import numpy as np
import math
import matplotlib.cm as cm
import cPickle as pic
import json
from collections import Counter
from glob import glob

#importing ddg pickle files
ddg_uqcon = pic.load(open("rosetta/out_pickles/ddg_dic_UQCON.pkl", "rb"))
ddg_cue = pic.load(open("rosetta/out_pickles/ddg_dic_CUE.pkl", "rb"))
ddg_monomer = pic.load(open("rosetta/out_pickles/ddg_dic_Monomer.pkl", "rb"))
ddg_otu = pic.load(open("rosetta/out_pickles/ddg_dic_OTU.pkl", "rb"))
ddg_rpn13 = pic.load(open("rosetta/out_pickles/ddg_dic_RPN13.pkl", "rb"))
ddg_sh3 = pic.load(open("rosetta/out_pickles/ddg_dic_SH3.pkl", "rb"))

#importing rosetta files
rosetta_test = pic.load(open("msd/msd_pickles/run_16_patch_5_propdic.pkl", "rb"))

#ddg of patch = sum of(ddg * fraction the mutation at that location shows up in rosetta)
#input (ddg dictionary, rosetta dictionary)
def ddg_calc(ddg_dict, rosetta_dict):
	ddg_list = []
	for aa_rosetta in rosetta_dict:
		total_ddg = 0
		for i in range (0,20):	
			try:
				total_ddg += (ddg_dict[int(aa_rosetta),(rosetta_dict[aa_rosetta][i][0])] * rosetta_dict[aa_rosetta][i][1])
			except KeyError:
				pass
		ddg_list.append(total_ddg)		
	return ddg_list

def entropy(s):
	p, lns = Counter(s), float(len(s))
	return -sum( count/lns * math.log(count/lns, 2) for count in p.values())

def entropy_ls(rosetta_dict):
	entropy_list = []
	for aa in rosetta_dict:
		aa_list = []
		for i in range(0, 20):
			aa_list.append(rosetta_dict[aa][i][1])
		entropy_list.append(entropy(aa_list))
	return entropy_list

def plotmaker(ddg_dict, rosetta_dict, name_string):
	x_val = ddg_calc(ddg_monomer, rosetta_test)
	y_val = entropy_ls(rosetta_test)
	plt.plot(x_val, y_val, "ro")
	plt.savefig('ddg_vs_entropy/' + name_string + '.pdf')
# print rosetta_test
# x_val = ddg_calc(ddg_monomer, rosetta_test)
# print len(rosetta_test)
# y_val = entropy_ls(rosetta_test)
# plt.plot(x_val, y_val, "ro")
#plt.savefig('somerthing.pdf')
# plt.show()
# print len(entropy_ls(rosetta_test))

filelist = glob('msd/msd_pickles/run*')
print filelist
for name in filelist:
	rosetta_file = pic.load(open(name, "rb"))
	plotmaker(ddg_uqcon, rosetta_file, name)

