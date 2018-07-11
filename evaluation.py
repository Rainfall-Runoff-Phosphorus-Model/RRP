# -*- coding: utf-8 -*-
"""
Rainfall-Runoff-Phosphorus model
(c) 2018, Agroscope, Sebastian Stoll

Evaluation of calibrated parameters + file generation for P model
"""

#%% Load libraries
from __future__ import division
import pandas as pd
from linecache import getline
import numpy as np
import matplotlib.pyplot as plt

#%% Generation of overland files
ovl=1

#%%use drainage information
dr=1

#%% Initialize Time
#Lag-factor
tau=1;

#Start at timestep...
st=0; #days
st=st*24; #hours

#%% Path of input & output files
pathI="...\\input\\"
pathO="...\\results\\"
pathOv="...\\overland\\"

#%% Load Calibration results
rNSE=np.loadtxt(pathO+'Calib_NSE.txt',delimiter=';') #Nash-Sutcliffe efficiencies of calibration
rP=np.loadtxt(pathO+'Calib_P.txt',delimiter=';') #Calibrated parameter sets

#%% Load Climate and discharge data
pName=pathI+'climate_evl.xlsx'
qName=pathI+'discharge_evl.xlsx'
clim=pd.read_excel(pName)
precip=clim['P(mm)'].values
evapot=clim['ETP(mm)'].values
qObs=pd.read_excel(qName)
qObsMM=qObs['Q(mm)'].values

#Datetime
dates=pd.to_datetime(clim.Year*10000+clim.Month*100+clim.Day,format='%Y%m%d')
dates=dates[st:len(precip)+1]


#%% load HRU and topo data

#read Header
hdr=[getline(pathI+'twi.txt',i) for i in range (1,7)]
values=[float(h.split(" ")[-1].strip()) \
    for h in hdr]
cols,rows,lx,ly,cell,nd=values
xres=cell
yres=cell*-1

#read topographic index
topo=np.loadtxt(pathI+'twi.txt',skiprows=6)
hru=np.loadtxt(pathI+'hru.txt',skiprows=6)
topo[topo==-9999]='NaN'
hru[hru==-9999]='NaN'

if dr==1:
    drain=np.loadtxt(pathI+'drain.txt',skiprows=6)
    drain[drain==-9999]='NaN'

#lamWell
lam1=topo[hru==3]
lam1[np.isnan(lam1)]=np.nanmean(lam1)
if dr==1:
    drain1=drain[hru==3]+0.5 #set TWI increase upon drainage
    drain1[drain1<1]=1    
    drain1[np.isnan(drain1)]=1
        
#lamPoor
lam2=topo[hru==4]
lam2[np.isnan(lam2)]=np.nanmean(lam2)
if dr==1:
    drain2=drain[hru==4]+0.5 #set TWI increase upon drainage
    drain2[drain2<1]=1 
    drain2[np.isnan(drain2)]=1
    
#lamUrban
lam3=topo[hru==2]
lam3[np.isnan(lam3)]=np.nanmean(lam1)
#lamForest
lam4=topo[hru==1]

l1=len(lam1)
l2=len(lam2)
l3=len(lam3)
l4=len(lam4)

max1=max(lam1)
max2=max(lam2)
min1=min(lam1)
min2=min(lam2)

tot=l1+l2+l3+l4


#%% Parameter vector
P=rP
aNSE=rNSE*0
sizeP=np.shape(P)

#%% Looping over parameter array
for k in xrange(0,sizeP[1]):
#for k in xrange(0,10):


    #%% Initialize States
    S1=np.zeros(1)
    S2=np.zeros(1)
    S4=np.zeros(1)
    SM1=np.zeros(1)
    SM2=np.zeros(1)
    SM4=np.zeros(1)
    dS1=np.zeros(1)
    dS2=np.zeros(1)
    dS4=np.zeros(1)
    qS1=np.zeros(len(precip)-st)
    qS2=np.zeros(len(precip)-st)
    qS4=np.zeros(len(precip)-st)
    qF1=np.zeros(len(precip)-st)
    qF2=np.zeros(len(precip)-st)
    qF3=np.zeros(len(precip)-st)
    lam01=np.zeros(1)
    lam02=np.zeros(1)
    A1F=np.zeros(1)
    A2F=np.zeros(1)
    Q=np.zeros((len(precip)-st,1))
    
    if ovl==1:
        fHRU1=np.int_(np.zeros((len(precip)-st,len(lam1))))
        fHRU2=np.int_(np.zeros((len(precip)-st,len(lam2))))
    
    
    #%% First time Step
    si=0.1 #initial soil water content
            
    # Soil water storage
    S1[0,]=P[0,k]*si #Initial SWS
    S2[0,]=P[1,k]*si #Initial SWS
    S4[0,]=P[12,k]*si #Initial SWS
    
    SM1[0,]=S1/P[0,k]; #Soil moisture
    SM2[0,]=S2/P[1,k]; #Soil moisture
    SM4[0,]=S4/P[12,k]; #Soil moisture
    
    # Slow flow HRU1 & 2 & 4
    qS1[0,]=SM1[0,]*P[8,k]; #slow flow
    qS2[0,]=SM2[0,]*P[9,k]; #slow flow
    qS4[0,]=SM4[0,]*P[13,k]; #slow flow
    
    # lamda threshold
    lam01[0,]=min1+((max1-min1)*(1-SM1)**P[10,k]);
    lam02[0,]=min2+((max2-min2)*(1-SM2)**P[11,k]);
    
    # Lamda areas
    if dr==1:
        A1F=np.sum(lam1*drain1>lam01)/l1
        A2F=np.sum(lam2*drain2>lam02)/l2
    else:    
        A1F=np.sum(lam1>lam01)/l1
        A2F=np.sum(lam2>lam02)/l2
       
    if ovl==1:
        fHRU1[0,lam1>lam01]=1
        fHRU2[0,lam2>lam02]=1
    
   
    # Fast flow HRU1 & 2
    qF1[0]=P[5,k]*precip[st-tau]*A1F
    qF2[0]=P[6,k]*precip[st-tau]*A2F
    
    # Fast flow HRU3 (urban)
    qF3[0]=P[7,k]*precip[st-tau]
    
    #deltaS
    dS1[0,]=precip[st]-qS1[0]-qF1[0]-evapot[st]
    dS2[0,]=precip[st]-qS2[0]-qF2[0]-evapot[st]
    dS4[0,]=precip[st]-qS4[0]-evapot[st]
    
    Q[0,]=(qS1[0]+qF1[0])*l1/tot + (qS2[0]+qF2[0])*l2/tot + qF3[0]*l3/tot + qS4[0]*l4/tot
    

    
    #%% Looping over time
    for d in xrange(1,len(precip)-st):
             
        #Soil moisture              
        S1[0,]=S1+dS1;
        S1[S1<0]=0;
        S2[0,]=S2+dS2
        S2[S2<0]=0
        S4[0,]=S4+dS4
        S4[S4<0]=0
           
        SM1[0,]=S1/P[0,k]
        SM1[SM1>1]=1;
        SM2[0,]=S2/P[1,k]
        SM2[SM2>1]=1;
        SM4[0,]=S4/P[12,k]
        SM4[SM4>1]=1;
        
        # Slow flow HRU1 & 2 & 4
        qS1[d]=SM1[0,]*P[8,k]; #slow flow
        qS2[d]=SM2[0,]*P[9,k]; #slow flow
        qS4[d]=SM4[0,]*P[13,k]; #slow flow
        
        # lamda threshold
        lam01[0,]=min1+((max1-min1)*(1-SM1)**P[10,k]);
        lam02[0,]=min2+((max2-min2)*(1-SM2)**P[11,k]);
        
        # Lamda areas
        A1F=np.sum(lam1>lam01)/l1
        A2F=np.sum(lam2>lam02)/l2
        
        if ovl==1:
            fHRU1[d,lam1>lam01]=1
            fHRU2[d,lam2>lam02]=1
        
        # Fast flow HRU1 & 2
        qF1[d]=P[2,k]*qF1[d-1] + P[5,k]*precip[st+d-tau]*A1F
        qF2[d]=P[3,k]*qF2[d-1] + P[6,k]*precip[st+d-tau]*A2F
        
        # Fast flow HRU3 (urban)
        qF3[d]=P[4,k]*qF3[d-1] + P[7,k]*precip[st+d-tau]
        
        #deltaS
        dS1[0,]=precip[st+d]-qS1[d]-qF1[d]-evapot[st+d]
        dS2[0,]=precip[st+d]-qS2[d]-qF2[d]-evapot[st+d]
        dS4[0,]=precip[st+d]-qS4[d]-evapot[st+d]
        
        Q[d,]=(qS1[d]+qF1[d])*l1/tot + (qS2[d]+qF2[d])*l2/tot + qF3[d]*l3/tot + qS4[d]*l4/tot
        

        
    #%% Aggregate hourly results to daily                  
    resultsQ=pd.DataFrame(Q)
    resultsQ=resultsQ.set_index(dates)
   
    if k==0:
        dailyQ=resultsQ.resample('d',how='sum') 
    else:
        dailyQ[k]=resultsQ.resample('d',how='sum')  
    dx=dailyQ.index
    
    
    
    #%% Calculate Nash-Sutcliffe
    Qobs=qObsMM[st/24:len(qObsMM)]
    meanQ=np.mean(Qobs)
   
    nse=1-(np.sum((Qobs-np.transpose(dailyQ[k].values))**2)/np.sum((Qobs-meanQ)**2))
    aNSE[k]=nse
    
       
    #%%Save results
    if ovl==1:
        np.savez_compressed(pathOv+'ovl_'+str(k),fHRU1=fHRU1,fHRU2=fHRU2,qS1=qS1,qS2=qS2,qF1=qF1,qF2=qF2,qF3=qF3,Q=Q)
  



   
#%% Plotting
    
dQ=dailyQ.values    
maxU=np.max(dQ,axis=1) 
minU=np.min(dQ,axis=1)
meanU=np.mean(dQ,axis=1)
dx=qObs.Date[int(st/24):len(qObsMM)]
di=dx.index.values
x=pd.date_range(start=dx[di[0]],end=dx[di[-1]],freq="D")

resultsP=pd.DataFrame(precip[st:len(precip)])
resultsP=resultsP.set_index(dates)
dailyP=resultsP.resample('d',how='sum')
dailyP.columns=['P']
width=1/1.5


plt.figure(1)
plt.subplot(211)
plt.bar(x,dailyP.P.values)
plt.ylabel('Daily P(mm)')

plt.subplot(212)
p1,=plt.plot(x,Qobs,'ro-')
plt.ylabel('Daily Discharge (mm)')
p2=plt.fill_between(x,minU,maxU,alpha=0.5,edgecolor='#1B2ACC',facecolor='#089FFF')
p3,=plt.plot(x,meanU,'k',color='#1B2ACC')
plt.legend([p1,p3],["observed","modelled"])


NSE=[rNSE,aNSE]
fig = plt.figure(2)
ax=fig.add_subplot(111)
bp=ax.boxplot(NSE)
ax.set_xticklabels(['Calibration','Evaluation'])
plt.ylabel('NSE')


print("Evaluation successful!")
print("Mean Nash-Sutcliffe Calibration: "+str(np.mean(rNSE)))
print("Mean Nash-Sutcliffe Evaluation: "+str(np.mean(aNSE)))
    
    
    
    
   
    
    
    
    
    
























        