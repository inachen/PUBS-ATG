import cPickle as pic


ubiquitin = "MQIFVKTLTGKTITLEVESSDTIDNVKSKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
newseq = "MQIFVKTLTGKTITLEVESSDTIDNVKSKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
wt_aa_dict = dict(((i, ubiquitin[i-1]),([],[])) for i in range (1,77))
wt_aa_dict[(77, 'STOP')]=([],[])

testset = dict(((i, newseq[i-1]),([float(i),float(i+1),float(i+2)],[float(i),float(i+1),float(i+2)])) for i in range(1,77))
testset[(77, 'STOP')]=([77.0,78.0,79.0],[77.0,78.0,79.0])

for i in range (1,76):
	testset[(i, newseq[i])] = ([i,i+1,i+2],[i,i+1,i+2])

count = 0

for key in wt_aa_dict.keys():
	if key in testset.keys():
		wt_aa_dict[key] = testset[key]
		count +=1
#dictionary now populated with wild-type sequences

wt_number_dict = dict((i, wt_aa_dict[i, ubiquitin[i-1]]) for i in range (1,77))
wt_number_dict[77] = ([77.0,77.0,77.0],[77.0,77.0,77.0])
#dictionary of wt seq with just position identifiers

pic.dump(wt_aa_dict, open("wt_aa_dict.pkl", "wb"))
pic.dump(wt_number_dict, open("wt_number_dict.pkl", "wb"))
pic.dump(testset, open("test_set.pkl", "wb"))