#!/usr/bin/env python3
# -*- coding:utf-8 -*-
###
# File: gamma.py
# Created Date: 2021-12-16 09:11:42
# Author: Nicolas zhan
# Contact: 1059291645@qq.com
# -----
# Last Modified: 2021-12-17 09:24:03
###
# %%
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import optimize as op
import os
import matplotlib.font_manager as font_manager
font_prop = font_manager.FontProperties(size=16)

# %%
def powerlaw(x,A,Gamma):
    return A*pow(x,-Gamma)

# %%
folder='to_disk_f'
try:
     os.mkdir(folder)
except FileExistsError:
    print("folder exits")
path=str('./'+folder+'/')

# %%
ebin=60
rbin=20
tbin=6400

e=np.genfromtxt("powerlog.dat")
dp=np.genfromtxt("Nrefl_inci.dat")
'''
e=np.genfromtxt("./e600/powerlog.dat")
dp=np.genfromtxt("./e600/refl_inci.dat")
'''
fp=np.zeros((ebin,rbin,tbin),dtype='double')
fp2=np.zeros((ebin,rbin,tbin),dtype='double')

# %%

idx=0
for k in range(tbin):
    rindx=0
    for j in range(rbin):
        for i in range(ebin):
            fp[i][j][k]=dp[k][2*i+1+rindx]
            idx+=1
        rindx+=2*ebin

# %%



# %% [markdown]
# fp [energy] [r_position] [time]

# %%
total=np.zeros(rbin)
for i in range(rbin):    
    for k in range(tbin):
        for j in range(ebin):
            total[i]+=fp[j][i][k]
            
rindex=np.linspace(50,500,20)
rindex=rindex[::-1]#序数 0 为最外环

fig,ax=plt.subplots()
ax.plot(rindex,total)
#plt.xlim(0.5,0.6)
plt.xlabel('radius')
plt.ylabel('counts')
plt.title('counts respect to radius')
plt.savefig(path+'counts_r.jpg')

# %%
total=np.zeros(tbin)
for k in range(tbin):
    for i in range(rbin):
        for j in range(ebin):
            total[k]+=fp[j][i][k]
            
rindex=np.linspace(50,500,20)
rindex=rindex[::-1]
tindx=np.linspace(0,100,tbin)
fig,ax=plt.subplots()
ax.plot(tindx,total)
#plt.xlim(0.5,0.6)
plt.xlabel('time')
plt.ylabel('counts')
plt.title('counts respect to time')
plt.savefig(path+'counts_t.jpg')

total=np.zeros((5,tbin))
for i in range(5):
    for k in range(tbin):
        for j in range(ebin):
            total[i,k]+=fp[j][3*i][k]
            
rindex=np.linspace(50,500,20)
rindex=rindex[::-1]
tindx=np.linspace(0,100,tbin)
fig,ax=plt.subplots()
for i in range(5):
    ax.plot(tindx,total[i,:],label='r'+str(3*i))
#plt.xlim(0.5,0.6)
plt.xlabel('time')
plt.ylabel('counts')
plt.legend(loc='bottom right')
plt.title('counts respect to t-diff-r')
plt.tight_layout()
plt.savefig(path+'count_diff-r.jpg')
