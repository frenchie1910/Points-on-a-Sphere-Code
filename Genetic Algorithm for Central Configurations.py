import numpy as np
import matplotlib.pyplot as plt
import random
from mpl_toolkits.mplot3d import Axes3D
import timeit
import math
from numpy.linalg import eig
import copy

start = timeit.default_timer()

#set population size
P=10

#set the number of points in each population member
n=22

# set the number of loops to run initial gradient flow for
loops=50

#set the number of generations you want to run the genetic algorithm for
gen=5

# set how often to output the energy

often=10

# Set initial time particle moves for
time=10

#open file to output energy during minimization 
filename="energygenetic"+str(n).zfill(3)
out = open(filename,'w')
out.close()

def makeNative(numpy):
    return getattr(numpy, "tolist", lambda: value)()

#the linear attractive force
def atforce(i):
    return -x[i]

#the inverse square repulsive force
def repforce(i):
    total=np.array([0,0,0])
    for j in range(n):
        if j==i:
            continue
        direction=x[i]-x[j]
        length=np.linalg.norm(direction)
        total=np.add(total,(direction/(length)**3))
    return total

#rotation matrix about z=0 plane if points have too many points
def Rz(theta):
  return np.matrix([[ np.cos(theta), -np.sin(theta), 0 ],
                   [ np.sin(theta), np.cos(theta) , 0 ],
                   [ 0           , 0            , 1 ]])

#creating an empty list to add the current population to
pop=[]

#creating a list to store the most recent population
popstore=[]


def POPatforce(a,i):
    return -np.array(pop[a][i])

def POPrepforce(a,i):
    total=np.array([0,0,0])
    for j in range(0,n):
        if j==i:
            continue
#         print(i,j)
#         print('pop',pop)
        #print('pop[a][i]',pop[a][i])
        #print('pop[a][j]',pop[a][j])
        direction=np.subtract(pop[a][i],pop[a][j])
        length=np.linalg.norm(direction)
        total=np.add(total,(direction/(length)**3))
    return total

#gradient flow function to run on each member of the population that returns the energy and coordinates in an array
def gradientflow(a):
    time=10
    for loop in range(0,30):
        old=copy.deepcopy(pop[a])
        for i in range(0,n):
            pop[a][i]=pop[a][i]+time*(POPrepforce(a,i)+POPatforce(a,i))
        # calculate the difference in energy
        #print(pop[a][i],a,g)
        difference=0.0
        energyrep=0.0
        energyat=0.0
        energyoldrep=0.0
        energyoldat=0.0
        for j in range(0,n):
            energyat+=(np.linalg.norm(pop[a][j])**2)/2
            energyoldat+=(np.linalg.norm(old[j])**2)/2
            for k in range(j+1,n):
                distance=np.sqrt(sum((pop[a][j]-pop[a][k])**2))
                energyrep+=1.0/distance
                distanceold=np.sqrt(sum((np.subtract(old[j],old[k]))**2))
                energyoldrep+=1.0/distanceold
        energy=energyat+energyrep
#         print('ENERGY',energy)
        energyold=energyoldat+energyoldrep
#         print('ENERGYOLD',energyold)
        difference=energy-energyold
        # accept or reject the move
        #print('difference',difference)
        if(difference<0.0):
            energy=energy
        else:
            pop[a]=old
            energy=energyold
            time=time*0.8
    return [[energy],[pop[a]]]

#loop to generate parent population
for p in range(0,P):
    random.seed()
    x=np.random.random((n,3))
#     print(x,p)
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
    pop.append(makeNative(x))

results=[]
#main loop that runs genetic algorithm
for g in range(0,gen):
    upper=[]
    lower=[]
    nextGEN=[]
    #seperate points above and below the z=0 plane for the parent configurations
    for j in range(0,P):
        xupper=[]
        xlower=[]
        for i in range(0,n):
            if pop[j][i][2]>0:
                xupper.append(pop[j][i])
            else:
                xlower.append(pop[j][i])
        upper.append(xupper)
        lower.append(xlower)
#[0] gets rid of 'array' and [0][i] gives the ith position in that list for the first element, [1][0] gives the first coordinate of the second parent
    #randomly mate them together
    for p in range(0,P):
        nextgen=[]
        ulen=len(upper)
        llen=len(lower)
        i=random.randint(0,ulen-1)
        j=random.randint(0,llen-1)
        nextgen=upper[i]+lower[j]
        #check how many points in each config
        if len(nextgen)!=n:
            if len(nextgen)<n:
                k=n-len(nextgen)
                for i in range(0,k):
                    nextgen.append(makeNative(np.random.random((3))))
            else:
                k=len(nextgen)-n
                random.shuffle(nextgen)
                for i in range(0,k):
                    nextgen.pop(i)
        nextGEN.append(nextgen)
        popstore=[]
        popstore.append([g,nextGEN]) #stores the last generation before the current one
    pop=nextGEN
#adds the current population to pop to then run the genetic algorithmon again
        #run gradient flow to relax configurations then run whole process again
#     pop=np.array(pop).ravel()
#     print('pop ravel',pop)
    for a in range(0,p):
        gradientflow(a)
        results.append(gradientflow(a))
finalenergy=[]
for i in range(p*g):
    finalenergy.append(results[i][0])

en=min(finalenergy)
print(round(en[0],4))

filename2= "energygenetic"+str(n).zfill(3)       
out1=open(filename2,'w')
out1 = open(filename,'a')
out1.write("{0:.4f} ".format(en[0]))
out1.close()

ind=np.argmin(finalenergy)
stop = timeit.default_timer()
execution_time = stop - start

print("Program Executed in "+str(execution_time))
#convert for 1 parent
# x1=[]
# x2=[]
# x3=[]
# for i in range(0,n):
#     x1.append(results[ind][i][0])
#     x2.append(results[ind][i][1])
#     x3.append(results[ind][i][2])
# 
# 
# #output final energy to the screen and points to a file
# # filename2= "pointsgenetic"+str(n).zfill(3)       
# # out=open(filename2,'w')
# # for i in range(0,n):
# #     for j in range(0,3):
# #         out.write("{0:.6f} ".format(pop[i,j]))
# #     out.write('\n')              
# # out.close()
#               
# 
# 
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(x1,x2,x3,color="purple",s=1000)
#            #s=np.pi*(4/3)1.65*3)
# ax.scatter(0,0,0,color='red',s=80)
# plt.axis('off')
# ax.view_init(azim=90.0, elev=90.0)
# 
# 
# stop = timeit.default_timer()
# execution_time = stop - start
# 
# print("Program Executed in "+str(execution_time))
# plt.show()