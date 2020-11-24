import sys
import numpy as np
import pandas as pd
#param = sys.stdin.readline().strip("\n").split(" ")
inf = sys.argv[1]
outf = sys.argv[2]

err = np.loadtxt(inf)
n = err.shape[0]
print(n)
rank = err.copy()
print(rank)
for i in range(n):
    tmp = err[i,1]
    #print(i,tmp)
    k = 0
    for j in range(n):
        if(j!=i and tmp > err[j,1]):
            k += 1
    #print(k)
    rank[i,1] = k
np.savetxt(outf ,rank, fmt="%.2d")
