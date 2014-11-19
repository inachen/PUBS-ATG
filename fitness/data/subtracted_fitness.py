import numpy as np
import sys

f1 = sys.argv[1]
f2 = sys.argv[2]

f1_mat = np.genfromtxt(sys.argv[1], delimiter=',')
f2_mat = np.genfromtxt(sys.argv[2], delimiter=',')

# HOW IT WURKS
#N,q = f1_mat.shape
#f1minf2 = np.zeros((N,q))
#
#for pos in range(N):
#    for am in range(q):
#        f1minf2[pos,am] = f1_mat[pos,am] - f2_mat[pos,am]
        
   
f_diff = f1_mat - f2_mat
np.savetxt("subtracted.csv", f_diff, delimiter=",")

