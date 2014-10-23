import cPickle as pic

ubiquitin = "MQIFVKTLTGKTITLEVESSDTIDNVKSKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
testset = pic.load(open("2500_filtered_seq.pkl", "rb"))
wildtypedict = dict(((i, ubiquitin[i]),([],[])) for i in range (1,76))
wildtypedict["(76, 'STOP')"]=([],[])

count = 0
for key in wildtypedict.keys():
	if key in testset.keys():
		wildtypedict[key] = testset[key]
		count +=1

print wildtypedict
print count