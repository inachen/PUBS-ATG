import sys
from Bio import SeqIO
import cPickle as pic
from barcode_differences.py import single_comp

barcode= str(raw_input("Barcode: ")).upper()
barcode = "CCGGTG"
filename = raw_input("sequencing fastq file name: ")
filename = "testseq.fastq"
newfile = open("CCGGTG.txt", "r+")
newpickle = pic.load(open("newpickle.pkl", "rb"))

count = 0

for record in SeqIO.parse(filename, "fastq"):
	"""pull out index from description"""
	indexread = record.description.split(":")[-1]
	indexreadcount = {indexread:0}
	"""record pulled out index for the barcode"""
	barcodecount = {barcode, indexreadcount}
	
	if single_comp(indexread, barcode) <3:
		"""count the number of reads pulled out"""
		count +=1
		#print record.description
		#newfile.write(record.description + str(record.letter_annotations["phred_quality"]))
		#print record.seq
		#newfile.write(record.seq)
		#print record.letter_annotations["phred_quality"]
		"""count the number of times the index is pulled out"""
		indexreadcount[indexread]+=1
		#newfile.write(record.letter_annotations["phred_quality"])
		#newfile.write(record.Seq)
		#newfile.write(record.phred)
		"""save the read to a new pickle"""
		pic.dump(record, open("newpickle.pkl", "wb"))
		#pic.dump(record.seq, open("newpickle.pkl", "wb"))
		#pic.dump(record.letter_annotations["phred_quality"], open("newpickle.pkl", "wb"))
print count
print barcodecount
newfile.close()
		

