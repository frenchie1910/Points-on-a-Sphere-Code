import numpy as np
import matplotlib.pyplot as plt
import random
from mpl_toolkits.mplot3d import Axes3D
import timeit
import math
from numpy.linalg import eig
import copy

#number of points
n=19

#create an empty array
DP=[]
for i in range(0,3):
    DP.append(0)
    
#open points file generated from gradient flow code
file='points'+str(n).zfill(3)
points=np.genfromtxt(open(file,'r'))

#calculate the dipole moment
def Dipole(n):
    for j in range(0,3):
        for i in range(0,n):
            DP[j]=DP[j]+points[i,j]
    return np.linalg.norm(DP)

print(Dipole(n))
print(DP)
filename="DipMom"+str(n).zfill(3)
out = open(filename,'w')
out.write('Dipole(n)')
out.close()