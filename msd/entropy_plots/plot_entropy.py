import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import collections
import csv
import numpy as np
import math
import matplotlib.cm as cm
from mpltools import style
from mpltools import color
import cPickle as pic
import pandas as pd

style.use('ggplot')

patch1_x = [0, 0.25, 0.5, 0.75, 1]
patch5_x = [0, 0.25, 0.5, 0.75]

names1 = ['M1', 'S19', 'S20', 'S57', 'N60', 'Q62', 'K63']
names5 = ['L8', 'R42', 'I44', 'H68', 'V70', 'R72']

patch1_1 = [0, 0, 0.279, 0.279, 0.279]
patch1_2 = [0, -0.444, -0.274, 0, 0]
patch1_3 = [0, -0.525, -0.525, -0.625, -0.625]
patch1_4 = [0, -0.217, -0.117, -0.217, -0.055]
patch1_5 = [0, -0.453, -0.823, -0.823, -0.823]
patch1_6 = [0, 0, -0.283, -0.561, -0.561]
patch1_7 = [0, -0.283, -0.283, -0.283, -0.283]

patch5_1 = [0, -1.234, -1.874, -1.799]
patch5_2 = [0, 0, 0, 0]
patch5_3 = [0, -0.392, -0.238, 0]
patch5_4 = [0, 1.169, 1.613, 1.613]
patch5_5 = [0, 0.138, 0.157, 0.320]
patch5_6 = [0, 0.374, 0.836, 0.836]

x_name = "Weight of other complexes"
y_name = "Change in Entropy"

fig = plt.figure(0)
plt.plot(patch1_x, patch1_1, marker = '.', linestyle='-', label = names1[0])
plt.plot(patch1_x, patch1_2, marker = '.', linestyle='-', label = names1[1])
plt.plot(patch1_x, patch1_3, marker = '.', linestyle='-', label = names1[2])
plt.plot(patch1_x, patch1_4, marker = '.', linestyle='-', label = names1[3])
plt.plot(patch1_x, patch1_5, marker = '.', linestyle='-', label = names1[4])
plt.plot(patch1_x, patch1_6, marker = '.', linestyle='-', label = names1[5])
plt.plot(patch1_x, patch1_7, marker = '.', linestyle='-', label = names1[6])
plt.xlabel(x_name)
plt.ylabel(y_name)
plt.title = "Patch 1 Entropy Plot"
lgd = plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
fig.savefig("patch1entropy.png", bbox_extra_artists=(lgd,), bbox_inches='tight')


fig = plt.figure(1)
plt.plot(patch5_x, patch5_1, marker = '.', linestyle='-', label = names5[0])
plt.plot(patch5_x, patch5_2, marker = '.', linestyle='-', label = names5[1])
plt.plot(patch5_x, patch5_3, marker = '.', linestyle='-', label = names5[2])
plt.plot(patch5_x, patch5_4, marker = '.', linestyle='-', label = names5[3])
plt.plot(patch5_x, patch5_5, marker = '.', linestyle='-', label = names5[4])
plt.plot(patch5_x, patch5_6, marker = '.', linestyle='-', label = names5[5])
plt.xlabel(x_name)
plt.ylabel(y_name)
plt.title = "Patch 5 Entropy Plot"
lgd = plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
fig.savefig("patch5entropy.png", bbox_extra_artists=(lgd,), bbox_inches='tight')



