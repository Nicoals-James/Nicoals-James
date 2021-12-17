#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from stingray import Lightcurve,AveragedCrossspectrum
import matplotlib.pyplot as plt
import math
# In[10]:


ne=60
nangle=3
nphi=25600
dt = 1/256  # seconds
exposure = 100 # seconds
#seg=exposure*256

enspec=0
psd=0
lag=1


fp2=np.zeros((ne,nangle,nphi))
dp2=np.genfromtxt('disk2.dat')
fp4=np.zeros((ne,nangle,nphi))
dp4=np.genfromtxt('comp2.dat')
fp6=np.zeros((ne,nangle,nphi))
dp6=np.genfromtxt('nesc2.dat')
fp8=np.zeros((ne,nangle,nphi))
dp8=np.genfromtxt('refl2.dat')

fp1=np.zeros((ne,nangle,nphi))
dp1=np.genfromtxt('disk1.dat')
fp3=np.zeros((ne,nangle,nphi))
dp3=np.genfromtxt('comp1.dat')
fp5=np.zeros((ne,nangle,nphi))
dp5=np.genfromtxt('nesc1.dat')
fp7=np.zeros((ne,nangle,nphi))
dp7=np.genfromtxt('refl1.dat')
for i in range(ne):
    id=0
    for j in range(nangle):
        for k in range(nphi):
            fp1[i][j][k]=dp1[i,3*id+1]#与输出结构有关
            fp2[i][j][k]=dp2[i,3*id+1]
            fp3[i][j][k]=dp3[i,3*id+1]
            fp4[i][j][k]=dp4[i,3*id+1]
            fp5[i][j][k]=dp5[i,3*id+1]#与输出结构有关
            fp6[i][j][k]=dp6[i,3*id+1]
            fp7[i][j][k]=dp7[i,3*id+1]
            fp8[i][j][k]=dp8[i,3*id+1]
            
            id+=1
Ntotal11=np.zeros(nphi)
Ntotal12=np.zeros(nphi)
Ntotal21=np.zeros(nphi)
Ntotal22=np.zeros(nphi)

E_s=1e-3
E_c=1e3
s_start = int(math.log10(E_s)*10+30)
s_end = int(math.log10(E_c)*10+30)
viewer=2#0是cosi [0~0.1]为平角，1为大约60°，2是正交观察cosi [0.9~1]

for i in range(nphi):
    for j in range(s_start,s_end):
        Ntotal11[i] +=fp3[j][viewer][i]#+fp5[j][viewer][i]+fp7[j][viewer][i]
        Ntotal12[i] +=fp4[j][viewer][i]#+fp6[j][viewer][i]+fp8[j][viewer][i]


E_s=1e-3
E_c=1e3
s_start = int(math.log10(E_s)*10+30)
s_end = int(math.log10(E_c)*10+30)
viewer=2#0是cosi [0~0.1]为平角，1为大约60°，2是正交观察cosi [0.9~1]

for i in range(nphi):
    for j in range(s_start,s_end):
        Ntotal21[i] +=fp7[j][viewer][i]
        Ntotal22[i] +=fp8[j][viewer][i]

        

for i in range(nphi):
    if Ntotal12[i] !=0:
        Ntotal11[i] = Ntotal11[i]/Ntotal12[i]
    if Ntotal22[i] !=0:
        Ntotal21[i] = Ntotal21[i]/Ntotal22[i]

times1 = np.arange(0, exposure, dt)  

lc1=Lightcurve(times1,Ntotal11)
lc2=Lightcurve(times1,Ntotal21)

avg_cs=AveragedCrossspectrum(lc2,lc1,5)
freq_lags, freq_lags_err = avg_cs.time_lag()

fig, ax = plt.subplots(1,1,figsize=(8,5))
ax.hlines(0, avg_cs.freq[0], avg_cs.freq[-1], color='black', linestyle='dashed', lw=2)
ax.errorbar(avg_cs.freq, freq_lags, yerr=freq_lags_err,fmt="o", lw=1, color='blue')
ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("Time lag (s)")
ax.tick_params(axis='x', labelsize=14)
ax.tick_params(axis='y', labelsize=14)
ax.tick_params(which='major', width=1.5, length=7)
ax.tick_params(which='minor', width=1.5, length=4)
plt.xscale('log')

plt.title("lag frequency dependence")
for axis in ['top', 'bottom', 'left', 'right']:
    ax.spines[axis].set_linewidth(1.5)
plt.savefig("cr-lag_whole-seg5.jpg")




