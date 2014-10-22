import sys

i=0
f = open("/data/ClassData-2014/data-oct16/lane1_Undetermined_L001_R1_001.fastq", "r")
copy = open("/data/atg/bruktrial/bigtestdata.fastq", "wt")

for i in range(1,10001):
	line = f.readline()
	copy.write(str(line))
f.close()
copy.close()


