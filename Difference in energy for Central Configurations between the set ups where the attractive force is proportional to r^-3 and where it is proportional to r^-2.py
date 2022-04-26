import numpy as np

results=[]
for n in range(2,21):
    file1='energycentral'+str(n).zfill(3)
    en1=np.genfromtxt(open(file1,'r'))
    file2='energythree'+str(n).zfill(3)
    en2=np.genfromtxt(open(file2,'r'))
    difference=en2-en1
    print(n,round(difference,4))
    results.append([n,round(difference,4)])

print(results)