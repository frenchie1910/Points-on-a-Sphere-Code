import numpy as np
import matplotlib.pyplot as plt
import random
from mpl_toolkits.mplot3d import Axes3D
import timeit
import math
from numpy.linalg import eig
import copy

start = timeit.default_timer()

# project one or all points to the sphere
def proj(x,j):
    if(j==n):
        for i in range(0,n):
            norm=np.sqrt(sum(x[i]**2))
            x[i]=x[i]/norm
    else:
        norm=np.sqrt(sum(x[j]**2))
        x[j]=x[j]/norm
    return x

#set the number of points
n=32

# set the number of loops
loops=50000000

# set how often to output the energy 
often=10
# maximum time particle moves for
time=10

#set the intitial temperature
temp=0.1

#set the rate of temperature reduction
reduc=0.95

# #open file to output energy during minimization 
# filename="energy"+str(n).zfill(3)
# out = open(filename,'w')
# out.close()

looplist=[]
energylist=[]
forcelist=[]

#assign random start points on the sphere
random.seed()
x=proj((2.0*np.random.random((n,3))-1.0),n)

def force(i):
    total=np.array([0,0,0])
    for j in range(n):
        if j==i:
            continue
        direction=x[i]-x[j]
        length=np.linalg.norm(direction)
        total=np.add(total,(direction/(length)**3))
    return total
# the main loop to reduce the energy        

for loop in range(0,loops+1):
    old=copy.deepcopy(x)
    for i in range(0,n):
        x[i]=x[i]+time*force(i)
        x=proj(x,i)
    # calculate the difference in energy
    difference=0.0
    energy=0.0
    energyold=0.0
    for j in range(0,n):
        for k in range(j+1,n):
            distance=np.sqrt(sum((x[j]-x[k])**2))
            energy+=1.0/distance
            distanceold=np.sqrt(sum((old[j]-old[k])**2))
            energyold+=1.0/distanceold
    difference=energy-energyold
    # accept or reject the move 
    if(difference<0.0):
        energy=energy
    else:
        P=(np.exp(difference/temp))/(1+np.exp(difference/temp))
        randomnumber=np.random.random(1)
        if randomnumber<P:
            energy=energy
        else:
            x=old
            time=time*0.8     
    if(loop%often==0):
        temp=reduc*temp
        print("{0} {1:.6f}".format(loop,energy))
#         out = open(filename,'a')
#         out.write("{0} {1:.6f} \n".format(loop,energy))
#         out.close()
        looplist.append(loop)
        energylist.append(energy)
    if math.isclose(difference,0,abs_tol=0.00000001):
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
    
#rotate coordinates

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


#to order eigenvalues so rotation works
def ordeig(v,w):
    if math.isclose(w[0]-w[2],w[1]-w[2],abs_tol=0.01):
        return v
    if math.isclose(w[0], w[2],abs_tol = 0.01):
        w[[1,2]] = w[[2,1]]
        v[1,2] = v[2,1]
        return v
    if math.isclose(w[1], w[2],abs_tol = 0.01):
        w[[0,2]] = w[[2,0]]
        v[0,2] = v[2,0]
        return v
    else:
        return v

v=v.T
R=np.array([[0,0,-1],[0,1,0],[1,0,0]])
x1r,x2r,x3r=np.matmul(v,[x1,x2,x3])
x1n,x2n,x3n=np.matmul(R,[x1r,x2r,x3r])
x=np.array([x1n,x2n,x3n])

#output final energy to the screen and points to a file
# filename2= "points"+str(n).zfill(3)       
# points=open(filename2,'w')
# for i in range(0,n):
#     for j in range(0,3):
#         points.write("{0:.6f} ".format(x.T[i,j]))
#     points.write('\n')              
# points.close()
    
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(xs, ys, zs,  rstride=4, cstride=4, color='yellow', alpha=0.5, linewidth=0)
ax.scatter(x1n,x2n,x3n,color="blue",s=80)
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
# plt.show()