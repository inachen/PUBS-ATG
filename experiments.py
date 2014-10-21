from Bio import SeqIO as seqio

for record in seqio.parse(testseq.fastq, "fastq"):
	print record
	

