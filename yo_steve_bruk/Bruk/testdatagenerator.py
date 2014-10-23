import cPickle as pic

ubiquitin = "MQIFVKTLTGKTITLEVESSDTIDNVKSKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
newseq = "MQIFVKTLTGKTITLEVESSDTIDNVKSKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
wt_aa_dict = dict(((i, ubiquitin[i-1]),([float(0), float(0), float(0)],[float(0), float(0), float(0)])) for i in range (1,77))
wt_aa_dict[(77, 'STOP')]=([],[])
testset = dict(((i, newseq[i-1]),([float(i),float(i),float(i)],[float(i),float(i),float(i)])) for i in range(1,77))
testset[(77, 'STOP')]=([77.0,77.0,77.0],[77.0,77.0,77.0])

for i in range (1,76):
	testset[(i, newseq[i])] = ([i,i,i],[i,i,i])

count = 0
print "printing empty dictionary"
print wt_aa_dict

for key in wt_aa_dict.keys():
	if key in testset.keys():
		wt_aa_dict[key] = testset[key]
		count +=1
#dictionary now populated with wild-type sequences

print "printing filled (pos, aa) dictionary"
print wt_aa_dict

wt_number_dict = dict((i, wt_aa_dict[i, ubiquitin[i-1]]) for i in range (1,77))
wt_number_dict[77] = ([77.0,77.0,77.0],[77.0,77.0,77.0])
#dictionary of wt seq with just position identifiers

print "printing filled pos dictionary"
print wt_number_dict

pic.dump(wt_aa_dict, open("wt_aa_dict.pkl", "wb"))
pic.dump(wt_number_dict, open("wt_number_dict.pkl", "wb"))
pic.dump(testset, open("test_set.pkl", "wb"))