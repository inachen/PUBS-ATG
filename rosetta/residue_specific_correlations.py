import numpy as np
import cPickle as pic
#from scipy.stats.stats import pearsonr 
fitness_dict = pic.load(open("input/D1S3fitness_scores.pkl", "rb"))
ddg = pic.load(open("out_pickles/ddg_dic.pkl", "rb"))


hydro = ['M', 'C', 'I', 'L', 'Y', 'F', 'W']
polar = ['Q', 'P', 'N', 'A', 'T', 'S', 'V', 'G']
charged = ['E', 'D', 'K', 'R', 'H']
compiled = {}

### adds ddG values to fitness scores, i.e {(#, aa):(fitness score, ddG)}###
count = 0
for key in ddg:
	try:
		fitness_score = fitness_dict[key]
		ddG_score = ddg[key]
		compiled[key] = (fitness_score, ddG_score)
	except KeyError as error:
		count +=1
		#print error
		#print ddG_dict[key]
		#print key in fitness_dict.keys()
		pass
print "unmatched key count: " + str(count)	

### parses by type of amino acid###
aa_type = int(raw_input("input # of choice >> 1: hydrophobic, 2: polar, 3: charged, 4: fitness scores, 5: ddG ???"))

if aa_type == 1:
	bin_type= hydro
elif aa_type ==2:
	bin_type = polar
elif aa_type ==3:
	bin_type = charged
elif aa_type==4:
	score_type = 0
	cut_off1, cut_off2 = float(raw_input("begining of fitness range:")), float(raw_input("end of fitness range:"))
else:
	score_type = 1
	cut_off1, cut_off2 = float(raw_input("begining of ddG range:")), float(raw_input("end of ddG range:"))

### structure is {(#, aa):([list of fitnesses], [list of corresponding ddGs])}###
pairwise = dict((i,([],[])) for i in range (1,77))

if aa_type<4:
	for key in compiled:
		aa_number = key[0]
		aa_identity = key[1]
		if aa_identity in bin_type:
			pairwise[aa_number][0].append(compiled[key][0])
			pairwise[aa_number][1].append(compiled[key][1])

else:
	for key in compiled:
		aa_number = key[0]
		aa_identity = key[1]
		if cut_off1<compiled[key][score_type]<cut_off2:
			pairwise[aa_number][0].append(compiled[key][0])
			pairwise[aa_number][1].append(compiled[key][1])

print "printing pairwise lists"
print pairwise

### computes the correlation coefficients for fitness vs. ddG for each ubq residue, {#:corr_coefficient}
correlation = dict((i,()) for i in range (1,77))

for key in pairwise:
	a = np.array(pairwise[key][0])
	b = np.array(pairwise[key][1])
	try:
		correlation[key]= np.corrcoef(a,b)[1,0]
	except:
		pass

final_data = []

for key in correlation:
	final_data.append(correlation[key])

print "printing correlations"
print final_data

filename = raw_input("save file as: ")

with open(filename, 'w') as filename:
	for item in final_data:
		try:
			filename.write("%s\t" % item)
		except:
			pass











	
	