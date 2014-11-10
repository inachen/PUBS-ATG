import cPickle as pic
from enrichment_value import enrich_val

ubiquitin = "MQIFVKTLTGKTITLEVESSDTIDNVKSKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
newseq = "MQIFVKTLTGKTITLEVESSDTIDNVKSKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
wt_aa_dict = dict(((i, ubiquitin[i-1]),([],[])) for i in range (1,77))
wt_aa_dict[(77, 'STOP')]=([],[])
testset= pic.load(open("filtered_seq_perfect.pkl", "rb"))

for key in wt_aa_dict.keys():
	if key in testset.keys():
		wt_aa_dict[key] = testset[key]

	else:
		wt_aa_dict[key] = ([1.0,1.0,1.0],[1.0,1.0,1.0])
		print "no reads for: " + str(key)
#dictionary now populated with wild-type sequences

wt_number_dict = dict((i, wt_aa_dict[i, ubiquitin[i-1]]) for i in range (1,77))
wt_number_dict[77] = (wt_aa_dict[77, 'STOP'])
#dictionary of wt seq with just position identifiers

count = 0

for key in wt_number_dict:
	for i in range (0,2):
		for j in range (0,3):
			count += wt_number_dict[key][i][j]
			
print count


