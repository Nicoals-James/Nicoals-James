#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from stingray import Lightcurve,AveragedCrossspectrum
import matplotlib.pyplot as plt
import shutil
import os
# In[10]:

save_folder='counts'

ne=60
nangle=3
nphi=25600
dt = 1/256  # seconds
exposure = 100 # seconds
#seg=exposure*256
fp1=np.zeros((ne,nangle,nphi))
dp1=np.genfromtxt('disk2.dat')
fp2=np.zeros((ne,nangle,nphi))
dp2=np.genfromtxt('comp2.dat')
fp3=np.zeros((ne,nangle,nphi))
dp3=np.genfromtxt('nesc2.dat')
fp4=np.zeros((ne,nangle,nphi))
dp4=np.genfromtxt('refl2.dat')
for i in range(ne):
    id=0
    for j in range(nangle):
        for k in range(nphi):
            fp1[i][j][k]=dp1[i,3*id+1]#与输出结构有关
            fp2[i][j][k]=dp2[i,3*id+1]
            fp3[i][j][k]=dp3[i,3*id+1]
            fp4[i][j][k]=dp4[i,3*id+1]
            
            id+=1
Ntotal1=np.zeros(nphi)
Ntotal2=np.zeros(nphi)
Ntotal3=np.zeros(nphi)
Ntotal4=np.zeros(nphi)
s_start=0 
s_end=60
viewer=2#0是cosi [0~0.1]为平角，1为大约60°，2是正交观察cosi [0.9~1]

for i in range(nphi):
    for j in range(s_start,s_end):
        Ntotal1[i] +=fp1[j][viewer][i]#按时间分，光子数
        Ntotal2[i] +=fp2[j][viewer][i]
        Ntotal3[i] +=fp3[j][viewer][i]
        Ntotal4[i] +=fp4[j][viewer][i]


path=str('./'+save_folder)
try:
    os.mkdir(path)
except FileExistsError:
    print("该文件夹已经存在，无需再创建")     
path=str(path+'/')

times1 = np.arange(0, exposure, dt)  
fig,ax=plt.subplots(1,1)
ax.plot(times1,Ntotal1,label='disk')
ax.legend(loc='upper right')
plt.savefig(path+'d-real-counts.jpg')  

fig,ax=plt.subplots(1,1)
ax.plot(times1,Ntotal2,label='comp')
plt.savefig(path+'c-real-counts.jpg') 

fig,ax=plt.subplots(1,1)
ax.plot(times1,Ntotal3,label='nesc')
plt.savefig(path+'n-real-counts.jpg') 

fig,ax=plt.subplots(1,1)
ax.plot(times1,Ntotal4,label='refl')
plt.savefig(path+'r-real-counts.jpg') 




