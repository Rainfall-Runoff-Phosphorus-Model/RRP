# -*- coding: utf-8 -*-
"""
Rainfall-Runoff-Phosphorus model
(c) 2018, Agroscope, Sebastian Stoll

Calibration of hydrological parameters
"""


#%% Load libraries

from __future__ import division
import pandas as pd
import numpy as np
import time
import matplotlib.pyplot as plt

#%% Calibration specifications

#use drainage information
dr=1

#save results
sv=1

#Path of input & output files
pathI="...\\input\\"
pathO="...\\results\\"

#Monte-Carlo Runs per Iteration
nMC=5000;

#Number of Iterations
nI=4;

#Lag-factor
tau=1;

#Start at timestep...
st=0; #days
st=st*24; #hours

#behavioral threshold
tN=0.8;

#%% Load Climate and discharge data

pName=pathI+'climate_cal.xlsx'
qName=pathI+'discharge_cal.xlsx'
clim=pd.read_excel(pName)
precip=clim['P(mm)'].values
evapot=clim['ETP(mm)'].values
qObs=pd.read_excel(qName)
qObsMM=qObs['Q(mm)'].values

#Datetime
dates=pd.to_datetime(clim.Year*10000+clim.Month*100+clim.Day,format='%Y%m%d')
dates=dates[st:len(precip)+1]

#%% load HRU and topo data

#read topographicI/HRU/drain
topo=np.loadtxt(pathI+'twi.txt',skiprows=6)
hru=np.loadtxt(pathI+'hru.txt',skiprows=6)
     
topo[topo==-9999]='NaN'
hru[hru==-9999]='NaN'

if dr==1:
    drain=np.loadtxt(pathI+'drain.txt',skiprows=6)
    drain[drain==-9999]='NaN'
    
#TWI well drained soils
lam1=topo[hru==3]
lam1[np.isnan(lam1)]=np.nanmean(lam1)
if dr==1:
    drain1=drain[hru==3]+0.5 #set TWI increase upon drainage
    drain1[drain1<1]=1    
    drain1[np.isnan(drain1)]=1
        
#WI poorly drained soils
lam2=topo[hru==4]
lam2[np.isnan(lam2)]=np.nanmean(lam2)
if dr==1:
    drain2=drain[hru==4]+0.5 #set TWI increase upon drainage
    drain2[drain2<1]=1 
    drain2[np.isnan(drain2)]=1
    
#TWI Urban
lam3=topo[hru==2]
lam3[np.isnan(lam3)]=np.nanmean(lam1)

#TWI Forest
lam4=topo[hru==1]

#length and maximum values of TWI
l1=len(lam1)
l2=len(lam2)
l3=len(lam3)
l4=len(lam4)

max1=max(lam1)
max2=max(lam2)
min1=min(lam1)
min2=min(lam2)

tot=l1+l2+l3+l4


#%% Loop over MC-Iteration

for k in xrange(0,nI):
    
    
    #%% Set parameters for MC-Iteration
    P=np.zeros((14,nMC))
    # SMAX: Max. SWS (mm) in HRU1 (well drained soil) 
    P[0,]=(800-100)*np.transpose(np.random.rand(nMC,1))+100
    # SMAX: Max. SWS (mm) in HRU2 (poorly drained soil) 
    P[1,]=(800-100)*np.transpose(np.random.rand(nMC,1))+100
    # SMAX: Max. SWS (mm) in HRU4 (forest) min=50, max=500
    P[12,]=(800-100)*np.transpose(np.random.rand(nMC,1))+100
    
        # a: Fast flopw decline rate in HRU1 (well drained soil) 
    P[2,]=(1-0)*np.transpose(np.random.rand(nMC,1))
    # a: Fast flopw decline rate in HRU2 (poorly drained soil) 
    P[3,]=(1-0)*np.transpose(np.random.rand(nMC,1))
    # a: Fast flopw decline rate in HRU3 (urban area) 
    P[4,]=(1-0)*np.transpose(np.random.rand(nMC,1))
    
    # b: Prop. of RF inti fast flow in HRU1 (well drained soil) 
    P[5,]=(0.7-0)*np.transpose(np.random.rand(nMC,1))
    # b: Prop. of RF inti fast flow in HRU2 (poorly drained soil)
    P[6,]=(1-0.3)*np.transpose(np.random.rand(nMC,1))+0.3
    # b: Prop. of RF inti fast flow in HRU3 (urban area) 
    P[7,]=(1-0.5)*np.transpose(np.random.rand(nMC,1))+0.5
    
    # c: Flow rate between SWS and slow flow comp. in HRU1 (well drained soil) 
    P[8,]=(1-0.3)*np.transpose(np.random.rand(nMC,1))+0.3
    # c: Flow rate between SWS and slow flow comp. in HRU2 (poorly drained soil) 
    P[9,]=(0.7-0)*np.transpose(np.random.rand(nMC,1))
    # c: Flow rate between SWS and slow flow comp. in HRU4 (forest) 
    P[13,]=(1-0.3)*np.transpose(np.random.rand(nMC,1))+0.3
    
    # n: Expansion control of areas cont. FF in HRU1 (well drained soil) 
    P[10,]=(10-1)*np.transpose(np.random.rand(nMC,1))+1
    # n: Expansion control of areas cont. FF in HRU2 (poorly drained soil) 
    P[11,]=(10-1)*np.transpose(np.random.rand(nMC,1))+1
            
        
    #%% Initialize States and variables
    S1=np.zeros((1,nMC))
    S2=np.zeros((1,nMC))
    S4=np.zeros((1,nMC))
    SM1=np.zeros((1,nMC))
    SM2=np.zeros((1,nMC))
    SM4=np.zeros((1,nMC))
    dS1=np.zeros((1,nMC))
    dS2=np.zeros((1,nMC))
    dS4=np.zeros((1,nMC))
    qS1=np.zeros((1,nMC))
    qS2=np.zeros((1,nMC))
    qS4=np.zeros((1,nMC))
    qF1=np.zeros((2,nMC))
    qF2=np.zeros((2,nMC))
    qF3=np.zeros((2,nMC))
    lam01=np.zeros((1,nMC))
    lam02=np.zeros((1,nMC))
    A1F=np.zeros((1,nMC))
    A2F=np.zeros((1,nMC))
    Q=np.zeros((len(precip)-st,nMC))

    
    
    #%%  Calculate 1. time step
    si=0.15 #initial soil water content
    
    # Soil water storage
    S1[0,]=P[0,]*si#Initial SWS
    S2[0,]=P[1,]*si #Initial SWS
    S4[0,]=P[12,]*si #Initial SWS
    
    SM1[0,]=S1/P[0,]; #Soil moisture
    SM2[0,]=S2/P[1,]; #Soil moisture
    SM4[0,]=S2/P[12,]; #Soil moisture
    
    # Slow flow HRU1 & 2 & 4
    qS1[0,]=SM1[0,]*P[8,]; #slow flow
    qS2[0,]=SM2[0,]*P[9,]; #slow flow
    qS4[0,]=SM4[0,]*P[13,]; #slow flow
    
    # lamda threshold
    lam01[0,]=min1+((max1-min1)*(1-SM1)**P[10,]);
    lam02[0,]=min2+((max2-min2)*(1-SM2)**P[11,]);
    
    # Lamda areas
    if dr==1:
       for l in xrange(0,nMC):
           A1F[0,l]=np.sum(lam1*drain1>lam01[0,l])/l1
           A2F[0,l]=np.sum(lam2*drain2>lam02[0,l])/l2
    else:    
        for l in xrange(0,nMC):
            A1F[0,l]=np.sum(lam1>lam01[0,l])/l1
            A2F[0,l]=np.sum(lam2>lam02[0,l])/l2
    
    # Fast flow HRU1 & 2
    qF1[1,]=P[5,]*precip[st-tau]*A1F
    qF2[1,]=P[6,]*precip[st-tau]*A2F
    
    # Fast flow HRU3 (urban)
    qF3[1,]=P[7,]*precip[st-tau]
    
    #deltaS
    dS1[0,]=precip[st]-qS1-qF1[1,]-evapot[st]
    dS2[0,]=precip[st]-qS2-qF2[1,]-evapot[st]
    dS4[0,]=precip[st]-qS4-evapot[st]
    
    Q[0,]=(qS1+qF1[1,])*l1/tot + (qS2+qF2[1,])*l2/tot + qF3[1,]*l3/tot + qS4*l4/tot
    
    qF1[0,]=qF1[1,]
    qF2[0,]=qF2[1,]
    qF3[0,]=qF3[1,]
    
    
    #%%  Calculate other time-steps
    for d in xrange(1,len(precip)-st):
             
        #Soil moisture
        S1[0,]=S1+dS1;
        S1[S1<0]=0;
        S2[0,]=S2+dS2
        S2[S2<0]=0
        S4[0,]=S4+dS4
        S4[S4<0]=0
           
        SM1[0,]=S1/P[0,]
        SM1[SM1>1]=1;
        SM2[0,]=S2/P[1,]
        SM2[SM2>1]=1;
        SM4[0,]=S4/P[12,]
        SM4[SM4>1]=1;
        
        # Slow flow HRU1 & 2
        qS1[0,]=SM1[0,]*P[8,]; #slow flow
        qS2[0,]=SM2[0,]*P[9,]; #slow flow
        qS4[0,]=SM4[0,]*P[13,]; #slow flow
    
        # lamda threshold
        lam01[0,]=min1+((max1-min1)*(1-SM1)**P[10,]);
        lam02[0,]=min2+((max2-min2)*(1-SM2)**P[11,]);
        
        
        # Lamda areas
        if dr==1:
           for l in xrange(0,nMC):
               A1F[0,l]=np.sum(lam1*drain1>lam01[0,l])/l1
               A2F[0,l]=np.sum(lam2*drain2>lam02[0,l])/l2
        else:    
            for l in xrange(0,nMC):
                A1F[0,l]=np.sum(lam1>lam01[0,l])/l1
                A2F[0,l]=np.sum(lam2>lam02[0,l])/l2
        
        # Fast flow HRU1 & 2
        qF1[1,]=P[2,]*qF1[0,] + P[5,]*precip[st+d-tau]*A1F
        qF2[1,]=P[3,]*qF2[0,] + P[6,]*precip[st+d-tau]*A2F
        
        # Fast flow HRU3 (urban)
        qF3[1,]=P[4,]*qF3[0,] + P[7,]*precip[st+d-tau]
        
        #deltaS
        dS1[0,]=precip[st+d]-qS1-qF1[1,]-evapot[st+d]
        dS2[0,]=precip[st+d]-qS2-qF2[1,]-evapot[st+d]
        dS4[0,]=precip[st+d]-qS4-evapot[st+d]
        
        Q[d,]=(qS1+qF1[1,])*l1/tot + (qS2+qF2[1,])*l2/tot + qF3[1,]*l3/tot
        
        qF1[0,]=qF1[1,]
        qF2[0,]=qF2[1,]
        qF3[0,]=qF3[1,]
        
     
                     
    resultsQ=pd.DataFrame(Q)
    resultsQ=resultsQ.set_index(dates)
    dailyQ=resultsQ.resample('d',how='sum')
    dailyQ=dailyQ.dropna(axis=0,how='all')
    
    #%%  Calculate Nash-Sutcliffe
    Qobs=qObsMM[st/24:len(qObsMM)]
    meanQ=np.mean(Qobs)
    Qobs=np.vstack([Qobs]*nMC).transpose()
    
    nse=1-(np.sum((Qobs-dailyQ.values)**2,axis=0)/np.sum((Qobs-meanQ)**2,axis=0))
    iS=np.where(nse>tN)[0]
    
  
    if k==0:
        rP=P[:,iS]
        rNSE=nse[iS]
    else:
        rP=np.concatenate((rP,P[:,iS]),axis=1)
        rNSE=np.concatenate((rNSE,nse[iS]))
    
    del Q,P,dailyQ,resultsQ




#%% Print and save calibration results



if sv==1:
    np.savetxt(pathO+'Calib_P.txt',rP,delimiter=';')
    np.savetxt(pathO+'Calib_NSE.txt',rNSE,delimiter=';')  

plt.figure(3)
ax = plt.subplot(531)
ax.set_title("SMAX 1")
plt.scatter(rP[0,],rNSE)
plt.ylabel('NSE')

ax=plt.subplot(532)
ax.set_title("SMAX 2")
plt.scatter(rP[1,],rNSE)
plt.ylabel('NSE')

ax=plt.subplot(533)
ax.set_title("SMAX 4")
plt.scatter(rP[12,],rNSE)
plt.ylabel('NSE')

ax=plt.subplot(534)
ax.set_title("FFLOW D 1")
plt.scatter(rP[2,],rNSE)
plt.ylabel('NSE')

ax=plt.subplot(535)
ax.set_title("FFLOW D 2")
plt.scatter(rP[3,],rNSE)
plt.ylabel('NSE')

ax=plt.subplot(536)
ax.set_title("FFLOW D 3")
plt.scatter(rP[4,],rNSE)
plt.ylabel('NSE')

ax=plt.subplot(537)
ax.set_title("PropFF 1")
plt.scatter(rP[5,],rNSE)
plt.ylabel('NSE')

ax=plt.subplot(538)
ax.set_title("PropFF 2")
plt.scatter(rP[6,],rNSE)
plt.ylabel('NSE')

ax=plt.subplot(539)
ax.set_title("PropFF 3")
plt.scatter(rP[7,],rNSE)
plt.ylabel('NSE')

ax=plt.subplot(5,3,10)
ax.set_title("PropSF 1")
plt.scatter(rP[8,],rNSE)
plt.ylabel('NSE')

ax=plt.subplot(5,3,11)
ax.set_title("PropSF 2")
plt.scatter(rP[9,],rNSE)
plt.ylabel('NSE')

ax=plt.subplot(5,3,12)
ax.set_title("PropSF 4")
plt.scatter(rP[13,],rNSE)
plt.ylabel('NSE')

ax=plt.subplot(5,3,13)
ax.set_title("EC FF1")
plt.scatter(rP[10,],rNSE)
plt.ylabel('NSE')

ax=plt.subplot(5,3,14)
ax.set_title("EC FF2")
plt.scatter(rP[11,],rNSE)
plt.ylabel('NSE')


print("Calibration done! Out of " + str(nMC*(nI+1)) + ", " +str( np.size(rP,1))+ " parameter combinations were accepted.")
        









