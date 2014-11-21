# plot heatmap from csv file

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

fpath = 'fitness_csv/Hydroxyurea_sub_DMSO_fitness.csv'

mat = np.loadtxt(open(fpath,"rb"),delimiter=",")

plt.clf()
fig = plt.figure()

heatmap = plt.pcolor(mat)

fig.savefig('Hydroxyurea_sub_DMSO.png')