#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from stingray import Lightcurve,AveragedCrossspectrum
import matplotlib.pyplot as plt

# In[10]:


ne=60
nangle=3
nphi=25600
dt = 1/256  # seconds
exposure = 100 # seconds
#seg=exposure*256
fp1=np.zeros((ne,nangle,nphi))
dp1=np.genfromtxt('comp1.dat')
fp2=np.zeros((ne,nangle,nphi))
dp2=np.genfromtxt('comp2.dat')
fp3=np.zeros((ne,nangle,nphi))
dp3=np.genfromtxt('disk1.dat')
fp4=np.zeros((ne,nangle,nphi))
dp4=np.genfromtxt('disk2.dat')
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

        


times1 = np.arange(0, exposure, dt)  
for i in range(nphi):
    if Ntotal2[i]!=0:
        Ntotal1[i]=Ntotal1[i]/Ntotal2[i]
        
for i in range(nphi):
    if Ntotal4[i]!=0:
        Ntotal3[i]=Ntotal3[i]/Ntotal4[i]
        
        
lc1=Lightcurve(times1,Ntotal1[:])       
lc2=Lightcurve(times1,Ntotal3[:])       
idx=0       
for i in range(nphi):
    if lc1.counts[i]==0 :
        #b=rnd.randint(2,10)/5
        #s2[i]=b
        idx+=1
print(idx)

fig,ax=plt.subplots(1,1)
ax.plot(times1,lc1.counts)
plt.savefig('comp-counts.jpg')

fig,ax=plt.subplots(1,1)
ax.plot(times1,lc2.counts)
plt.savefig('disk-counts.jpg')

# In[ ]:
''''''
avgps=AveragedCrossspectrum(lc2,lc1,5)
freq_lags, freq_lags_err = avgps.time_lag()

fig, ax = plt.subplots(1,1,figsize=(8,5))
ax.hlines(0, avgps.freq[0], avgps.freq[-1], color='black', linestyle='dashed', lw=2)
ax.hlines(0.0276,avgps.freq[0], avgps.freq[-1], color='black', linestyle='dashed', lw=2)
ax.errorbar(avgps.freq, freq_lags, yerr=freq_lags_err,fmt="o", lw=1, color='blue')
#ax.plot(avgps.freq,freq_lags)
ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("Time lag (s)")
ax.tick_params(axis='x', labelsize=14)
ax.tick_params(axis='y', labelsize=14)
ax.tick_params(which='major', width=1.5, length=7)
ax.tick_params(which='minor', width=1.5, length=4)
plt.xscale("log")
plt.xlim(0,50)
for axis in ['top', 'bottom', 'left', 'right']:
    ax.spines[axis].set_linewidth(1.5)
plt.savefig("dc-lag-sum-seg5.jpg")


avgps1=AveragedCrossspectrum(lc2,lc1,1)
freq_lags1, freq_lags_err1 = avgps1.time_lag()

fig, ax = plt.subplots(1,1,figsize=(8,5))
ax.hlines(0, avgps1.freq[0], avgps1.freq[-1], color='black', linestyle='dashed', lw=2)
ax.hlines(0.0276,avgps1.freq[0], avgps1.freq[-1], color='black', linestyle='dashed', lw=2)
ax.errorbar(avgps1.freq, freq_lags1, yerr=freq_lags_err1,fmt="o", lw=1, color='blue')
#ax.plot(avgps.freq,freq_lags)
ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("Time lag (s)")
ax.tick_params(axis='x', labelsize=14)
ax.tick_params(axis='y', labelsize=14)
ax.tick_params(which='major', width=1.5, length=7)
ax.tick_params(which='minor', width=1.5, length=4)
plt.xscale("log")
plt.xlim(0,50)
for axis in ['top', 'bottom', 'left', 'right']:
    ax.spines[axis].set_linewidth(1.5)
plt.savefig("dc-lag-sum-seg1.jpg")
# In[ ]:





# In[ ]:




