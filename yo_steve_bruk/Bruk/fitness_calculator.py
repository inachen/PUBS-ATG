import cPickle as pic
wt_aa_dict = pic.load(open("wt_aa_dict.pkl", "rb"))
wt_number_dict = pic.load(open("wt_number_dict.pkl", "rb"))
testset = pic.load(open("test_set.pkl", "rb"))

print testset[57, 'D']
print testset[57, 'S']
print wt_number_dict[57]

for key in testset:
	aa_pos = key[0]	
	for i in range (0,3):
		testset[key][0][i] /= wt_number_dict[aa_pos][0][i]
		testset[key][1][i] /= wt_number_dict[aa_pos][1][i]
# adds wt counts for aa position to mutant counts for aa position in testset

print testset[57, 'D']
print testset[57, 'S']

"""
ubiquitin = "MQIFVKTLTGKTITLEVESSDTIDNVKSKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"

wt_aa_dict = dict(((i, ubiquitin[i-1]),([],[])) for i in range (1,77))
wt_aa_dict[(77, 'STOP')]=([],[])

testset = open("picklename.pkl", "rb")

count = 0

for key in wt_aa_dict.keys():
	if key in testset.keys():
		wt_aa_dict[key] = testset[key]
		count +=1
		
#dictionary now populated with wild-type sequences

wt_number_dict = dict((i, wt_aa_dict[i, ubiquitin[i-1]]) for i in range (1,77))
wt_number_dict[77] = ([77.0,77.0,77.0],[77.0,77.0,77.0])

#dictionary of wt seq with just position identifiers
"""