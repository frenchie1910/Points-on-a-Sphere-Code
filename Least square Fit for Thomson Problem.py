import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

n=30

xdata=[]
fitydata=[]
for i in range(2,n+1):
    xdata.append(i)

ydata=[0.5,1.73205081,3.67423461,6.47469149,9.98528137,14.45297741,19.67528786,25.75998653,32.71694946,40.59645051,49.16525306,58.85323061,69.30636330,80.67024411,92.91165630,106.05040483,120.08446745,135.08946756,150.88156833,167.64162240,185.28753615,203.93019066,223.34707405,243.81276030,265.13332632,287.30261503,310.49154236,334.63443992,359.60394590]
def fit(n,b):
    return 0.5*n**2+b*n**(3/2)

c,cov=curve_fit(fit,xdata,ydata)

print(c)
for i in range(2,n+1):
        fitydata.append(fit(i,c[0]))

plt.xlabel("Number of points")
plt.ylabel("Energy")

plt.plot(xdata,ydata,color='lightcoral', linewidth=4, linestyle=':')
plt.plot(xdata,fitydata)

plt.show()

#plot the difference between actual energy value and function

plt.xlabel("Number of points")
plt.ylabel("Difference in Energy")

differencedata=np.subtract(fitydata,ydata)
plt.plot(xdata,differencedata)
plt.xticks(range(1, 31))

plt.show()
