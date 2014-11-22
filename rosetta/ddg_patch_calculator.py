import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import collections
from collections import namedtuple
import csv
import numpy as np
import math
import matplotlib.cm as cm
import cPickle as pic
import json

#importing ddg pickle files
ddg_uqcon = pic.load(open("rosetta/out_pickles/ddg_dic_UQCON.pkl", "rb"))
ddg_cue = pic.load(open("rosetta/out_pickles/ddg_dic_CUE.pkl", "rb"))
ddg_monomer = pic.load(open("rosetta/out_pickles/ddg_dic_Monomer.pkl", "rb"))
ddg_otu = pic.load(open("rosetta/out_pickles/ddg_dic_OTU.pkl", "rb"))
ddg_rpn13 = pic.load(open("rosetta/out_pickles/ddg_dic_RPN13.pkl", "rb"))
ddg_sh3 = pic.load(open("rosetta/out_pickles/ddg_dic_SH3.pkl", "rb"))

#importing rosetta files
rosetta_test = pic.load(open("msd/msd_pickles/run_1_patch_1_propdic.pkl", "rb"))

#ddg of patch = sum of(ddg * fraction the mutation at that location shows up in rosetta)
#input (ddg dictionary, rosetta dictionary)
def ddg_calc(ddg_dict, rosetta_dict):
	total_ddg = 0
	for aa_rosetta in rosetta_dict:
		for i in range (0,20):
			try:
				total_ddg += (ddg_dict[int(aa_rosetta),(rosetta_dict[aa_rosetta][i][0])] * rosetta_dict[aa_rosetta][i][1])
			except KeyError:
				pass
	return total_ddg

print rosetta_test
#print ddg_calc(ddg_monomer, rosetta_test)
# print ddg_sh3[int('30'), rosetta_test['30'][0][0]]