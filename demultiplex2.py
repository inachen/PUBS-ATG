import sys

workfile = open("newtestseq.fastq", "r")
newfile = open("ATCACG.fastq", "r+")

for line in workfile.readlines():
	barcodeline = str(line.split(":")[-1])
	barcode = barcodeline.split()[0]
	if barcode =="ATCACG":
		newfile.write(line)
		
	print barcode

	
workfile.close()
newfile.close()
		