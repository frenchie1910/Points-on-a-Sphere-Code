#plotting figures for r3

import numpy as np
import matplotlib.pyplot as plt
import random
from mpl_toolkits.mplot3d import Axes3D
import timeit
import math
from numpy.linalg import eig
import copy

for n in range(2,21):
    file='pointsthree'+str(n).zfill(3)
    points=np.genfromtxt(open(file,'r'))

    x1=[]
    x2=[]
    x3=[]
    x1s=[]
    x2s=[]
    x3s=[]

    data=points.T
    for i in range(0,n):
        x1.append(data[0,i])
        x2.append(data[1,i])
        x3.append(data[2,i])
        
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(0,0,0,color="red",s=200)
    ax.scatter(x1,x2,x3,color="purple",s=300)
    ax.scatter(x1s,x2s,x3s,color="green",s=1000)

               #s=np.pi*(4/3)*1.65**3)
    # ax.scatter(0,0,0,color='red',s=80)
    plt.axis('off')
    ax.view_init(azim=90.0, elev=90.0)

    # for i in range(0,n):
    #     for j in range(0,n):
    #         ax.plot([x1[i],x1[j]],[x2[i],x2[j]],[x3[i],x3[j]])
plt.show()
