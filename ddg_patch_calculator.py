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

patches = {'patch_1':['1','19','20','57','60','62','63'],
'patch_2':['42','44','45','46','48','49'],
'patch_3': ['6','8','42','44','66','68','70','73'],
'patch_4': ['44', '45', '57', '58', '59', '60', '62', '65'], 
'patch_5': ['8', '42', '44', '68', '70', '72']}

#ddg of patch = sum of(ddg * fraction the mutation at that location shows up in rosetta)
#input (ddg dictionary, rosetta dictionary)
def ddg_calc(ddg_dict, rosetta_dict, name_string):
	ddg_list = []
	for aa in patches[name_string]:
		total_ddg = 0
		for i in range (0,20):	
			try:
				total_ddg += (ddg_dict[int(aa),(rosetta_dict[aa][i][0])] * rosetta_dict[aa][i][1])
			except KeyError:
				pass
		ddg_list.append(total_ddg)		
	return ddg_list

def entropy(s):
	p, lns = Counter(s), float(len(s))
	return -sum( count/lns * math.log(count/lns, 2) for count in p.values())

def entropy_ls(rosetta_dict, name_string):
	entropy_list = []
	for aa in patches[name_string]:
		aa_list = []
		for i in range(0, 20):
			aa_list.append(rosetta_dict[aa][i][1])
		entropy_list.append(entropy(aa_list))
	return entropy_list

def plotmaker(ddg_dict, rosetta_dict, name_string):
	split_name = name_string.split('_')
	compose_name = split_name[3] +'_'+ split_name[4]
	x_val = ddg_calc(ddg_dict, rosetta_dict, compose_name)
	y_val = entropy_ls(rosetta_dict, compose_name)
	plt.plot(x_val, y_val, "ro")
	plt.xlabel('DDG Values')
	plt.ylabel('Shannon Entropy')
	plt.title('Correlation of DDG and Shannon Entropy')
	plt.savefig('ddg_vs_entropy/' + name_string[16:37] + '.pdf')
	print name_string
	print x_val
	print y_val
# print rosetta_test
# x_val = ddg_calc(ddg_monomer, rosetta_test)
# print len(rosetta_test)
# y_val = entropy_ls(rosetta_test)
# plt.plot(x_val, y_val, "ro")
#plt.savefig('somerthing.pdf')
# plt.show()
# print len(entropy_ls(rosetta_test))

#importing rosetta files
filelist = glob('msd/msd_pickles/run*')
for name in filelist:
	rosetta_file = pic.load(open(name, "rb"))
	plotmaker(ddg_uqcon, rosetta_file, name)

#file1 = pic.load(open('msd/msd_pickles/run_1_patch_1_propdic.pkl', "rb"))
#file1 = filelist[0]
