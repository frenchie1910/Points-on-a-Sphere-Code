import numpy as np
import math

symno=[]
n=19
file='pointscentralI'+str(n).zfill(3)
points=np.genfromtxt(open(file,'r'))

def rotmat(theta):
    c=np.cos(theta)
    s=np.sin(theta)
    Rz=np.array([[c,-s,0],[s,c,0],[0,0,1]])
    Ry=np.array([[c,0,s],[0,1,0],[-s,0,c]])
    Rx=np.array([[1,0,0],[0,c,-s],[0,s,c]])
    return Rz

#print('points',points)
def rotpoints(j):
    i=0
    k=0
    xrot=np.matmul(rotmat(2*math.pi/(12-j)),points.T)
    xrot=xrot.T
    while k<n:
        if math.isclose(np.sqrt((xrot[i,0]-points[k,0])**2+(xrot[i,1]-points[k,1])**2+(xrot[i,2]-points[k,2])**2),0,abs_tol=0.1):
            i=i+1
            k=0
            if i==n:
                print('the configuration has a C', 12-j, 'symmetry')
                exit()
        else:
            k=k+1

for m in range(0,11):
    rotpoints(m)