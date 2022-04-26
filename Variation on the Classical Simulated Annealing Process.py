import numpy as np
import matplotlib.pyplot as plt
import random
from mpl_toolkits.mplot3d import Axes3D
from numpy.linalg import eig
import math
import timeit

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
n=5

# set the number of loops
loops=1000000

# set how often to output the energy 
often=100

#set the intitial temperature
temp=0.1

#set the rate of temperature reduction
reduc=0.5

#set the max amplitude of the angle
k=100
 
#open file to output energy during minimization 
filename="energyannealing"+str(n).zfill(3)
out = open(filename,'w')
out.close()

looplist=[]
energylist=[]

#assign random start points on the sphere
random.seed()
x=proj((2.0*np.random.random((n,3))-1.0),n)

# calculate the initial energy
energy=0.0
for i in range(0,n):
    for j in range(i+1,n):
        distance=np.sqrt(sum((x[i]-x[j])**2))
        energy+=1.0/distance

def randrot(alpha, beta, gamma):
    a=np.array([np.cos(alpha)*np.cos(beta),np.cos(alpha)*np.sin(beta)*np.sin(gamma)-np.sin(alpha)*np.cos(gamma),np.cos(alpha)*np.sin(beta)*np.cos(gamma)+np.sin(alpha)*np.sin(gamma)])
    b=np.array([np.sin(alpha)*np.cos(beta),np.sin(alpha)*np.sin(beta)*np.sin(gamma)+np.cos(alpha)*np.cos(gamma),np.sin(alpha)*np.sin(beta)*np.cos(gamma)-np.cos(alpha)*np.sin(gamma)])
    c=np.array([-np.sin(beta),np.cos(beta)*np.sin(gamma),np.cos(beta)*np.cos(gamma)])
    MAT=np.array([a.flatten(),b.flatten(),c.flatten()])
    return MAT

# the main loop to reduce the energy        

start = timeit.default_timer()
for loop in range(0,loops+1):
    # randomly choose a point to move
    i=random.randint(0,n-1)
    alpha=(np.random.random(1))/k
    beta=(np.random.random(1))/k
    gamma=(np.random.random(1))/k
    randomnumber=np.random.random(1)

    # store the old coordinates of this point
    old=np.array(x[i])
    # randomly move this point
    x[i]=np.matmul(randrot(alpha, beta, gamma),x[i].T)
    x=proj(x,i)

    # calculate the difference in energy
    difference=0.0
    for j in range(0,n):
        if(j!=i):
            distance=np.sqrt(sum((x[i]-x[j])**2))
            distanceold=np.sqrt(sum((old-x[j])**2))
            difference=difference+1.0/distance-1.0/distanceold;

    # accept or reject the move 
    if(difference<0.0):
        energy=energy+difference
    else:
        P=(np.exp(difference/temp))/(1+np.exp(difference/temp))
        if randomnumber<P:
            energy=energy+difference
        else:
            x[i]=old

    # output energy to screen and a file
    
    if(loop%often==0):
        temp=reduc*temp
        print("{0} {1:.6f}".format(loop,energy))
        out = open(filename,'a')
        out.write("{0} {1:.6f} \n".format(loop,energy))
        out.close()
        looplist.append(loop)
        energylist.append(energy)
        k=k+1
    if temp<0.01:
        break
stop = timeit.default_timer()
execution_time = stop - start
# output final energy to the screen and points to a file
print("Final energy = {0:.6f} \n".format(min(energylist)))
print(execution_time)
filename2= "pointsannealing"+str(n).zfill(3)       
points=open(filename2,'w')
for i in range(0,n):
    for j in range(0,3):
        points.write("{0:.6f} ".format(x[i,j]))
    points.write('\n')              
points.close()

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


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(xs, ys, zs,  rstride=4, cstride=4, color='yellow', alpha=0.5, linewidth=0)
#ax.scatter(x1,x2,x3,color="black",s=80) #old coordinates before rotation
ax.scatter(x1n,x2n,x3n,color="blue",s=80)

# #drawing principal axis
# if math.isclose(w[0], w[1],abs_tol = 0.01):
#     if math.isclose(w[1], w[2],abs_tol = 0.01):
#         line = np.linspace(-1, 1, 2)
#         ax.plot3D(v[0,0]*line, v[0,1]*line, v[0,2]*line, 'gray')
#     else:
#         line = np.linspace(-1, 1, 2)
#         ax.plot3D(v[2,0]*line, v[2,1]*line, v[2,2]*line, 'green')
# if math.isclose(w[1], w[2],abs_tol = 0.01):
#     line = np.linspace(-1, 1, 2)
#     ax.plot3D(v[0,0]*line, v[0,1]*line, v[0,2]*line, 'blue')
    


# line = np.linspace(-1, 1, 2)
# ax.plot3D(v[0,0]*line, v[0,1]*line, v[0,2]*line, 'gray')
# line = np.linspace(-1, 1, 2)
# ax.plot3D(v[2,0]*line, v[2,1]*line, v[2,2]*line, 'green')
# line = np.linspace(-1, 1, 2)
# ax.plot3D(v[1,0]*line, v[1,1]*line, v[1,2]*line, 'purple')


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

# plt.show()