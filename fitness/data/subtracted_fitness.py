import numpy as np
import sys


f1_mat = np.genfromtxt(sys.argv[1], delimiter=',')
#HU data
f2_mat = np.genfromtxt(sys.argv[2], delimiter=',')
#other group data

# HOW IT WURKS
N,q = f1_mat.shape
f1minf2 = np.zeros((N,q))

for pos in range(N):
    for am in range(q):
        f1minf2[pos,am] = f1_mat[pos,am] - f2_mat[pos,am]
        if f1_mat[pos,am] == 0. or f2_mat[pos,am] == 0.:
            f1minf2[pos,am] == 0.
   
#f_diff = f1_mat - f2_mat
np.savetxt("{0}.csv".format(sys.argv[3]), f1minf2, delimiter=",")
#ouput csv file

import os
import matplotlib.pyplot as plt
plt.pcolormesh(f1minf2, cmap='Spectral')
plt.savefig("{0}.png".format(sys.argv[3]))
