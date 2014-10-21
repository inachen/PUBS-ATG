import bio
from bio import seqIO

record_iterator = seqIO.parse("testsmall.fastq", "fastq")
first_record = next(record_iterator)

print first_record
