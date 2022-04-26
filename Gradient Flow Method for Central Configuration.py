import numpy as np
import matplotlib.pyplot as plt
import random
from mpl_toolkits.mplot3d import Axes3D
import timeit
import math
from numpy.linalg import eig
import copy

start = timeit.default_timer()

#set the number of points
n=58

# set the number of loops
loops=100000

# set how often to output the energy 
often=100
# maximum time particle moves for
time=10

# #open file to output energy during minimization 
# filename="energycentral"+str(n).zfill(3)
# out = open(filename,'w')
# out.close()

looplist=[]
energylist=[]
forcelist=[]

#assign random start points on the sphere
random.seed()
x=2.0*np.random.random((n,3))-1.0

def atforce(i):
    return -x[i]

def repforce(i):
    total=np.array([0,0,0])
    for j in range(n):
        if j==i:
            continue
        direction=x[i]-x[j]
        length=np.linalg.norm(direction)
        total=np.add(total,(direction/(length)**3))
    return total

for loop in range(0,loops+1):
    old=copy.deepcopy(x)
    for i in range(0,n):
        x[i]=x[i]+time*(repforce(i)+atforce(i))
    # calculate the difference in energy
    difference=0.0
    energyrep=0.0
    energyat=0.0
    energyoldrep=0.0
    energyoldat=0.0
    for j in range(0,n):
        energyat+=(np.linalg.norm(x[j])**2)/2
        energyoldat+=(np.linalg.norm(old[j])**2)/2
        for k in range(j+1,n):
            distance=np.sqrt(sum((x[j]-x[k])**2))
            energyrep+=1.0/distance
            distanceold=np.sqrt(sum((old[j]-old[k])**2))
            energyoldrep+=1.0/distanceold
    energy=energyat+energyrep
    energyold=energyoldat+energyoldrep
    difference=energy-energyold
    # accept or reject the move 
    if(difference<0.0):
        energy=energy
    else:
        x=old
        energy=energyold
        time=time*0.8
    if(loop%often==0):
        print("{0} {1:.6f}".format(loop,energy))
        looplist.append(loop)
        energylist.append(energy)
    if math.isclose(difference,0,abs_tol=0.000001):
        break

print("{0} {1:.6f}".format(loop,energy))

#Create a sphere
theta, phi = np.mgrid[0.0:np.pi:100j, 0.0:2.0*np.pi:100j]
xs = np.sin(theta)*np.cos(phi)
ys = np.sin(theta)*np.sin(phi)
zs = np.cos(theta)

#convert data
x1=[]
x2=[]
x3=[]
for i in range(0,n):
    x1.append(x[i,0])
    x2.append(x[i,1])
    x3.append(x[i,2])

Ixx=[]
Iyy=[]
Izz=[]
Ixy=[]
Iyz=[]
Ixz=[]

def get_inertia_matrix(x1,x2,x3):
    # Moment of intertia tensor
    for i in range(0,n):
        Ixx.append(x2[i]**2 + x3[i]**2)
        Iyy.append(x1[i]**2 + x3[i]**2)
        Izz.append(x1[i]**2 + x2[i]**2)
        Ixy.append(x1[i] * x2[i])
        Iyz.append(x2[i] * x3[i])
        Ixz.append(x1[i] * x3[i])
    mIxx = np.sum(Ixx)
    mIyy = np.sum(Iyy)
    mIzz = np.sum(Izz)
    mIxy = -np.sum(Ixy)
    mIyz = -np.sum(Iyz)
    mIxz = -np.sum(Ixz)
    I = np.array([[mIxx, mIxy, mIxz], [mIxy, mIyy, mIyz], [mIxz, mIyz, mIzz]])
    return I

I = get_inertia_matrix(x1,x2,x3)

w,v=eig(I)
E=list(zip(w,v))
E.sort()


# #to order eigenvalues so rotation works
# def ordeig(v,w):
#     if math.isclose(w[0]-w[2],w[1]-w[2],abs_tol=0.01):
#         return v
#     if math.isclose(w[0], w[2],abs_tol = 0.01):
#         w[[1,2]] = w[[2,1]]
#         v[1,2] = v[2,1]
#         return v
#     if math.isclose(w[1], w[2],abs_tol = 0.01):
#         w[[0,2]] = w[[2,0]]
#         v[0,2] = v[2,0]
#         return v
#     else:
#         return v

v=v.T
xR=np.array([[1,0,0],[0,0,-1],[0,1,0]])
yR=np.array([[0,0,-1],[0,1,0],[1,0,0]])
zR=np.array([[0,-1,0],[1,0,0],[0,0,1]])
x1r,x2r,x3r=np.matmul(v,[x1,x2,x3])
x1n,x2n,x3n=np.matmul(yR.T,[x1r,x2r,x3r])
x=np.array([x1n,x2n,x3n])

#output final energy to the screen and points to a file
filename2= "pointscentralI"+str(n).zfill(3)       
out=open(filename2,'w')
for i in range(0,n):
    for j in range(0,3):
        out.write("{0:.6f} ".format(x.T[i,j]))
    out.write('\n')              
out.close()
# 
# filename2= "energy"+str(n).zfill(3)       
# out1=open(filename2,'w')
# out1 = open(filename,'a')
# out1.write("{0:.4f} ".format(energy))
# out1.close()              


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#ax.scatter(x1,x2,x3,color="blue",s=np.pi*(4/3)*1.65**3)
ax.scatter(x1n,x2n,x3n,color="blue",s=1000)
ax.scatter(0,0,0,color='red',s=80)
plt.axis('off')
ax.view_init(azim=90.0, elev=90.0)

# ax.set_xlim([-1.0,1.0])
# ax.set_ylim([-1.0,1.0])
# ax.set_zlim([-1.0,1.0])
# ax.set_aspect("auto")
# ax.set_title("{0} ".format(n)+"points on a sphere")
# 
# fig = plt.figure()
# ax2=fig.add_subplot(212)
# ax2.plot(looplist,energylist)
# ax2.set_xlabel("iteration number")
# ax2.set_ylabel("energy")

stop = timeit.default_timer()
execution_time = stop - start

print("Program Executed in "+str(execution_time))
plt.show()

