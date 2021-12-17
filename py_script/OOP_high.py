#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from stingray import Lightcurve, Crossspectrum, AveragedCrossspectrum
from scipy.interpolate import interp1d
import shutil
import os
import math
from stingray import Lightcurve, Powerspectrum, AveragedPowerspectrum
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
font_prop = font_manager.FontProperties(size=16)

class Load:
    def __init__(self,filename):
        self.filename = filename
    
    def light_read(self,ndat,ne,nangle,nphi):
        fp = np.zeros((ndat,ne,nangle,nphi))
        dp = []
        for i in range(ndat):
            dp.append(np.genfromtxt(self.filename[i]))
            for j in range(ne):
                id=0
                for k in range(nangle):
                    for l in range(nphi):
                        fp[i][j][k][l]=dp[-1][j,3*id+1]
                        id=id+1
        return fp

class Analysis(Load):
    def __init__(self,filename,ne=60,nangle=3,nphi=25600,dt =1/256):
        super().__init__(filename)
        exposure= nphi*dt
        self.times =  np.arange(0,exposure,dt)
        self.ne = ne
        self.nangle = nangle
        self.nphi = nphi
        self.ndat = len(filename)
        self.fp,self.er=self.light_read(self.ndat,self.ne,self.nangle,self.nphi)
    
    def light_read(self,ndat,ne,nangle,nphi):
        fp = np.zeros((ndat,ne,nangle,nphi))
        dp = []
        for i in range(ndat):
            dp.append(np.genfromtxt(self.filename[i]))
            for j in range(ne):
                id=0
                for k in range(nangle):
                    for l in range(nphi):
                        fp[i][j][k][l]=dp[-1][j,3*id+1]
                        id=id+1
        er = dp[-1][:,0]
        return fp,er
    
    def dir2save(self,path):
        self.path = str('/' + path)
        try:
            os.mkdir(path)
        except FileExistsError:
            print("该文件夹已经存在，无需再创建")
    
    def E_Spectrum(self,viewer):
        
        self.flux=np.zeros((self.ndat,self.ne))
        for i in range(self.ndat):
            for j in range(self.ne):
                for k in range(self.nphi):
                    self.flux[i][j] += self.fp[i][j][viewer][k]/self.nphi
        
        energy=self.er
        fig,ax=plt.subplots(1,1,figsize=(9,6))
        ax.loglog(energy,self.flux[0],c='blue',label='disk')
        ax.loglog(energy,self.flux[2],c='black',label='nesc')
        ax.loglog(energy,self.flux[4],c='red',label='comp')
        ax.loglog(energy,self.flux[6],c='purple',label='refl')
        ax.legend(loc='upper right')
        #plt.savefig('E_spectrum_high.jpg')

    def LightCurve(self,E_s,E_c,viewer):
        
        s_start = int(math.log10(E_s)*10+self.ne/2)
        s_end = int(math.log10(E_c)*10+self.ne/2)
        print(s_start,s_end)
        lgc=np.zeros((int(self.ndat/2),self.nphi))
        for i in range(int(self.ndat/2)):
            for j in range(s_start,s_end):
                for k in range(self.nphi):
                    try:
                        tmp =float(self.flux[2*i+1][j]*self.fp[i*2][j][viewer][k])/float(self.fp[2*i+1][j][viewer][k])
                    except ZeroDivisionError:
                        lgc[i][k] += 0
                    else:
                        lgc[i][k] += tmp

        return lgc

    def Powerspectrum(self,lgc,interval,tr,f_rebin):
        
        lc1 = Lightcurve(self.times, lgc[:])
        avg_ps1 = AveragedPowerspectrum(lc1,interval,norm='frac')
        avg_ps1 = avg_ps1.rebin_log(f=f_rebin)
        fig, ax1 = plt.subplots(1,1,figsize=(6,6))
        #ax1.plot(avg_ps1.freq,avg_ps1.power,c='purple',label='original')
        if tr:
            ax1.plot(avg_ps1.freq,avg_ps1.freq*avg_ps1.power,c='k')
        else:
            ax1.plot(avg_ps1.freq,avg_ps1.power,c='k',label='developed')
        #ax1.errorbar(avg_ps1.freq,avg_ps1.power,yerr=avg_ps1.power_err,c='green',ecolor='r',capsize=5,label='10')
        #ax1.errorbar(avg_ps2.freq,avg_ps2.power,yerr=avg_ps2.power_err,c='purple',ecolor='r',capsize=5,label='50')
        #ax1.errorbar(avg_ps3.freq,avg_ps3.power,yerr=avg_ps3.power_err,c='orange',ecolor='r',capsize=5,label='100')
        ax1.set_xlabel("Frequency (Hz)", fontproperties=font_prop)
        ax1.set_ylabel("Power (raw)", fontproperties=font_prop)
        ax1.set_yscale('log')
        ax1.tick_params(axis='x', labelsize=16)
        ax1.tick_params(axis='y', labelsize=16)
        ax1.tick_params(which='major', width=1.5, length=7)
        ax1.tick_params(which='minor', width=1.5, length=4)
        for axis in ['top', 'bottom', 'left', 'right']:
            ax1.spines[axis].set_linewidth(1.5)
    
        legend = plt.legend(loc='upper right', shadow=True, fontsize='x-large')
        plt.xscale('log')
        plt.title('avg_PSD')
        plt.savefig('avg_PSD_high.png')
    

    def LagFreq(self,lgc1,lgc2,interval,deltt):

        lc9 = Lightcurve(self.times, lgc1[:],dt=deltt)
        lc10 = Lightcurve(self.times,lgc2[:],dt=deltt)
        
        avg_cs = AveragedCrossspectrum(lc9, lc10,interval)
        '''
        fig, ax = plt.subplots(1,1,figsize=(10,6))
        ax.plot(lc9.time, lc9.counts, lw=2, color='blue')
        ax.plot(lc10.time,lc10.counts, lw=2, color='red')
        ax.set_xlabel("Time (s)", fontproperties=font_prop)
        ax.set_ylabel("Counts (cts)", fontproperties=font_prop)
        ax.tick_params(axis='x', labelsize=16)
        ax.tick_params(axis='y', labelsize=16)
        plt.title("lightcurve-comparison")
        plt.savefig("lightcurve-comparison.jpg")

        fig, ax = plt.subplots(1,1,figsize=(10,6))
        ax.plot(avg_cs.freq, avg_cs.power, lw=2, color='blue')
        
        plt.title("cross spectrum")
        plt.savefig("crossspectrum.jpg")
        '''
        freq_lags, freq_lags_err = avg_cs.time_lag()

        fig, ax = plt.subplots(1,1,figsize=(8,5))
        ax.hlines(0, avg_cs.freq[0], avg_cs.freq[-1], color='black', linestyle='dashed', lw=2)
        ax.errorbar(avg_cs.freq, freq_lags, yerr=freq_lags_err,fmt="o", lw=1, color='blue')
        ax.set_xlabel("Frequency (Hz)", fontproperties=font_prop)
        ax.set_ylabel("Time lag (s)", fontproperties=font_prop)
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)
        ax.tick_params(which='major', width=1.5, length=7)
        ax.tick_params(which='minor', width=1.5, length=4)
        plt.xscale('log')
        plt.title("lag frequency dependence")
        for axis in ['top', 'bottom', 'left', 'right']:
            ax.spines[axis].set_linewidth(1.5)
        plt.savefig("time-lag-high.jpg")


# In[2]:


file=['disk1.dat','disk2.dat','nesc1.dat','nesc2.dat','comp1.dat','comp2.dat','refl1.dat','refl2.dat']
exp1 = Analysis(file)


# In[4]:


#exp1.dir2save('figure')
exp1.E_Spectrum(2)
lgc1 = exp1.LightCurve(27,70,2)
lgc2 = exp1.LightCurve(70,150,2)
exp1.Powerspectrum(lgc1[0,:]+lgc1[1,:]+lgc1[2,:]+lgc1[3,:],100,1,0.02)
exp1.LagFreq(lgc1[0,:]+lgc1[1,:]+lgc1[2,:]+lgc1[3,:],lgc2[0,:]+lgc2[1,:]+lgc2[2,:]+lgc2[3,:],8,1/64)# if this value is positive - lgc2 lags lgc1









# In[ ]:




