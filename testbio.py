import sys
# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# import collections
# import csv
# import numpy as np
# import math
# import matplotlib.cm as cm
# from mpltools import style
# from mpltools import color
import cPickle as pic


from Bio import SeqIO

for record in SeqIO.parse(open("test.fastq", "rU"), "fastq"):
	print record.seq
	print record.description
	print record.letter_annotations['phred_quality']