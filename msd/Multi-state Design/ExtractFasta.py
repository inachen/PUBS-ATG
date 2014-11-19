from glob import glob
patchlist = glob("MSD*sh.o*")

for filename in patchlist:
	infile = open(filename, 'r')
	outfilename = filename + ".fasta"
	outfile = open(outfilename, 'w')
	count = 0
	for line in infile:
		if ": AA:" in line and count < 125:
			count +=1
			outstring = str(">%s_%d\n") % (filename, count)
			fields = line.split()
			for element in fields:
				if "AA:" in element:
					aa = element[-1]
					outstring = outstring + aa
			outstring = outstring +"\n"
			outfile.write(outstring)
	outfile.close()

