import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

n=30
    
xdata=[]
ydata=[]
fitydata=[]
for i in range(2,n+1):
    file='energycentral'+str(i).zfill(3)
    energycentral=np.genfromtxt(open(file,'r'))
    xdata.append(i)
    ydata.append(energycentral)
    print(energycentral)

    
def fit(n,a,b):
    return a*n**(5/3)-b*n

c,cov=curve_fit(fit,xdata,ydata)

print(c)
for i in range(2,n+1):
        fitydata.append(fit(i,c[0],c[1]))


plt.xlabel("Number of points")
plt.ylabel("Energy")

plt.plot(xdata,ydata,color='lightcoral', linewidth=4, linestyle=':')
plt.plot(xdata,fitydata)

plt.show()

#plot the difference between actual energy value and function


plt.xlabel("Number of points")
plt.ylabel("Energy")

differencedata=np.subtract(fitydata,ydata)
plt.plot(xdata,differencedata)

plt.show()