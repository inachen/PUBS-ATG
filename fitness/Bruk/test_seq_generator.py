import sys

i=0
j=1+4*int(raw_input("how many sequences: "))
f = open("/data/ClassData-2014/data-oct16/lane1_Undetermined_L001_R1_001.fastq", "r")
copy = open("/data/atg/bruktrial/testdata.fastq", "wt")



for i in range(1,j):
	line = f.readline()
	copy.write(str(line))
f.close()
copy.close()


