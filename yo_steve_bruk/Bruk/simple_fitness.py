import cPickle as pic
from enrichment_value import enrich_val

testset= pic.load(open("filtered_seq_perfect.pkl", "rb"))
wt_counts= pic.load(open("wt_counts.pkl", "rb"))
testset= enrich_val(testset)

for key in testset:
		for i in range (0,2):
			for j in range (0,3):
				testset[key][i][j] /= wt_counts[('enrich_wt')][i][j]
		
pic.dump(testset, open("fitness_scores.pkl", "wb"))	
