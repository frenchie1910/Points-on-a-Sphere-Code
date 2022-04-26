import numpy as np
import math

symno=[]
pointsnew=[]
n=5
file='points'+str(n).zfill(3)
points=np.genfromtxt(open(file,'r'))

# print(points)
for i in range(0,n):
    
    pointsnew[i,2]=-points[i,2]

# print('points ref',points)
    

#print('points',points)

i=0
k=0
if math.isclose(np.sqrt((xrot[i,0]-points[k,0])**2+(xrot[i,1]-points[k,1])**2+(xrot[i,2]-points[k,2])**2),0,abs_tol=0.1):
    i=i+1
    k=0
    if i==n:
        print('the configuration has a C', 12-j, 'symmetry')
        exit()
else:
    k=k+1