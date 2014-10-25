import cPickle as pic
import numpy as np

ubiquitin = "MQIFVKTLTGKTITLEVESSDTIDNVKSKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
wt_aa_dict = dict(((i, ubiquitin[i-1]),([],[])) for i in range (1,77))
wt_aa_dict[(77, 'STOP')]=([],[])
wt_aa_dict[(0, 'WT')]=([],[])
testset= pic.load(open("filtered_seq_perfect.pkl", "rb"))
num_dic = pic.load(open("aminotonumber.pkl", "rb"))


######################## for normalization with ALL wt counts #######################
#gets all the counts for WT sequences
for key in wt_aa_dict:
	if key in testset.keys():
		wt_aa_dict[key] = testset[key]

	else:
		wt_aa_dict[key] = ([0.0,0.0,0.0],[0.0,0.0,0.0])
		print "no reads for: " + str(key)
#dictionary now populated with wild-type sequences

#print wt_aa_dict

#####################################################################################

all_wt_counts = {'WT':([0.0,0.0,0.0],[0.0,0.0,0.0])}

#add up all wt counts
for key in wt_aa_dict:
	for i in range (0,2):
		for j in range (0,3):
			all_wt_counts['WT'][i][j] += wt_aa_dict[key][i][j]

#print all_wt_counts

#####################################################################################

total_counts = {'total':([0.0,0.0,0.0],[0.0,0.0,0.0])}

#add up all counts
for key in testset:
	for i in range (0,2):
		for j in range (0,3):
			total_counts['total'][i][j] += testset[key][i][j]

#print total_counts

#####################################################################################

norm_all_wt_counts = {'norm_wt':([0.0,0.0,0.0],[0.0,0.0,0.0])}

#normalize all_wt_count by total_count
for i in range (0,2):
	for j in range (0,3):
		norm_all_wt_counts['norm_wt'][i][j] = all_wt_counts['WT'][i][j]/total_counts['total'][i][j]

#print norm_all_wt_counts

enrich_all_wt_counts = {'enrich_wt':([0.0,0.0,0.0],[0.0,0.0,0.0])}

for i in range (0,2):
	for j in range (0,3):
		enrich_all_wt_counts['enrich_wt'][i][j] = norm_all_wt_counts['norm_wt'][i][j]/norm_all_wt_counts['norm_wt'][i][0]
			
#print enrich_all_wt_counts

pic.dump(enrich_all_wt_counts, open("all_wt_counts.pkl", "wb"))

####################### for normalization with 'WT' barcoded counts #######################

norm_wt_barcode_counts = {'norm_wt_barcode':([0.0,0.0,0.0],[0.0,0.0,0.0])}
#print norm_wt_barcode_counts

for i in range (0,2):
	for j in range (0,3):
		norm_wt_barcode_counts['norm_wt_barcode'][i][j]= testset[0, 'WT'][i][j]/total_counts['total'][i][j]

#print norm_wt_barcode_counts

enrich_wt_barcode_counts = {'enrich_wt':([0.0,0.0,0.0],[0.0,0.0,0.0])}

for i in range (0,2):
	for j in range (0,3):
		enrich_wt_barcode_counts['enrich_wt'][i][j]= norm_wt_barcode_counts['norm_wt_barcode'][i][j]/ norm_wt_barcode_counts['norm_wt_barcode'][i][0]
		

#print enrich_wt_barcode_counts

pic.dump(enrich_wt_barcode_counts, open("all_wt_barcode_counts.pkl", "wb"))


			
