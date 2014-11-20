import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import collections
from collections import namedtuple
import csv
import numpy as np
import math
import matplotlib.cm as cm
import copy
import os.path

from glob import glob
from Bio import SeqIO
import cPickle as pickle

dmso_in = open('Codon_FitScore_DMSO_day2.pkl', 'rb')
dmso_fit = pickle.load(infile)
dmso_in.close()

caffeine_in = open('Codon_FitScore_Caffeine_day2.pkl', 'rb')
caffeine_fit = pickle.load(infile)
caffeine_in.close()

hu_in = open('Codon_FitScore_HU_day2.pkl', 'rb')
hu_fit = pickle.load(infile)
hu_in.close()



print dmso_fit