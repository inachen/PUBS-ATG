import cPickle as pic
from enrichment_value import enrich_val

ubiquitin = "MQIFVKTLTGKTITLEVESSDTIDNVKSKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
wt_aa_dict = dict(((i, ubiquitin[i-1]),([],[])) for i in range (1,77))
wt_aa_dict[(77, 'STOP')]=([],[])
wt_aa_dict[(0, 'WT')]=([],[])
testset= pic.load(open("filtered_seq_perfect.pkl", "rb"))

for key in wt_aa_dict:
	if key in testset.keys():
		wt_aa_dict[key] = testset[key]

	else:
		wt_aa_dict[key] = ([0.0,0.0,0.0],[0.0,0.0,0.0])
		print "no reads for: " + str(key)
#dictionary now populated with wild-type sequences

print wt_aa_dict

wt_counts = {'WT':([0.0,0.0,0.0],[0.0,0.0,0.0])}

print wt_counts

#add up all wt counts
for key in wt_aa_dict:
	for i in range (0,2):
		for j in range (0,3):
			wt_counts['WT'][i][j] += wt_aa_dict[key][i][j]

print wt_counts

total_counts = {'total':([0.0,0.0,0.0],[0.0,0.0,0.0])}

#add up all counts
for key in testset:
	for i in range (0,2):
		for j in range (0,3):
			total_counts['total'][i][j] += testset[key][i][j]

print total_counts

norm_wt_counts = {'norm_wt':([0.0,0.0,0.0],[0.0,0.0,0.0])}

#normalize wt_count by total_count
for i in range (0,2):
	for j in range (0,3):
		norm_wt_counts['norm_wt'][i][j] = wt_counts['WT'][i][j]/total_counts['total'][i][j]

print norm_wt_counts

enrich_wt_counts = {'enrich_wt':([0.0,0.0,0.0],[0.0,0.0,0.0])}

for i in range (0,2):
	for j in range (0,3):
		enrich_wt_counts['enrich_wt'][i][j] = norm_wt_counts['norm_wt'][i][j]/norm_wt_counts['norm_wt'][i][0]
			
print enrich_wt_counts

testset= enrich_val(testset)





pic.dump(enrich_wt_counts, open("wt_counts.pkl", "wb"))
pic.dump(testset, open("test_set.pkl", "wb"))
