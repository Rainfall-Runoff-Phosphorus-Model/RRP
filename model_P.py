# -*- coding: utf-8 -*-
"""
Rainfall-Runoff-Phosphorus model
(c) 2018, Agroscope, Sebastian Stoll

Model Phosphorus Load/Concentration and generate Hydrological Risk Map
"""

#%% Load libraries
from __future__ import division
import pandas as pd
from linecache import getline
import numpy as np
import glob

#%% Path of input & output files
pathI="...\\input\\"
pathO="...\\results\\"
pathOv="...\\overland\\"

#%% load spatial input data

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

#lamWell
lam1=topo[hru==3]
lam1[np.isnan(lam1)]=np.nanmean(lam1)
       
#lamPoor
lam2=topo[hru==4]
lam2[np.isnan(lam2)]=np.nanmean(lam2)
    
#lamUrban
lam3=topo[hru==2]
lam3[np.isnan(lam3)]=np.nanmean(lam1)
#lamForest
lam4=topo[hru==1]

#soil phosphorus
sph=0 #soil phosphorus map available

if sph==1:
   sP=np.loadtxt(pathI+'sP.txt',skiprows=6)
   sP[sP==-9999]='NaN'
   pWell=sP[hru==3]
   pPoor=sP[hru==4]
else:
    sPc=18 # default water soluable P concentration (mg/kg)
    pWell=lam1*0+sPc
    pPoor=lam2*0+sPc
    
l1=len(lam1)
l2=len(lam2)
l3=len(lam3)
l4=len(lam4)

max1=max(lam1)
max2=max(lam2)
min1=min(lam1)
min2=min(lam2)

tot=l1+l2+l3+l4
area=(len(lam1)+len(lam2)+len(lam3)+len(lam4))*xres*xres


#%% Overland data
#number of overland files
nf=len(glob.glob(pathOv + '*'))

#Load first file
npzfile=np.load(pathOv+'ovl_0.npz')
qF1=npzfile['qF1']
 
#define P module parameters and initialize variables   
n=0.25 #new water fraction
DRPbase=0.05 #Baseflow concentration (mg/l)
t=5621 # Pick a certain timestep to create FastFlow map
LTot=np.zeros((len(qF1),nf))
QTot=np.zeros((len(qF1),nf))
h1=np.zeros((1,l1))
h2=np.zeros((1,l2))

#%% Q New and Old and corresponding P losses

for k in xrange(0,nf):
    
    npzfile=np.load(pathOv+'ovl_'+ str(k) + '.npz')    
    
    #size of overland variables
    s1=np.shape(npzfile['fHRU1'])
    s2=np.shape(npzfile['fHRU2'])
    
    #load overland variables
    fHRU1=npzfile['fHRU1']
    fHRU2=npzfile['fHRU2']
            
    h1=h1+fHRU1[t,:]
    h2=h2+fHRU2[t,:]
               
    #initialize variables           
    q1Old=np.zeros(np.shape(npzfile['fHRU1']))
    q2Old=np.zeros(np.shape(npzfile['fHRU2']))
    q1New=np.zeros(np.shape(npzfile['fHRU1']))
    q2New=np.zeros(np.shape(npzfile['fHRU2']))
    
    #read results from file    
    qF1=npzfile['qF1']
    qF2=npzfile['qF2']
    qS1=npzfile['qS1']
    qS2=npzfile['qS2']
    Q=npzfile['Q']
    
    del npzfile
    
    #% Areas with fast flow
    f1A=np.sum(fHRU1,1)/s1[1]
    f2A=np.sum(fHRU2,1)/s2[1]
    
    #%Old and new water
    for l in xrange(0,s1[0]):
        q1New[l,]=n*qF1[l]*(f1A[l]**(-1))*fHRU1[l,]
        q1Old[l,]=(1-n)*qF1[l]*(f1A[l]**(-1))*fHRU1[l,]+qS1[l]
    
    q1Old=np.nan_to_num(q1Old)
    q1New=np.nan_to_num(q1New)
    
    for l in xrange(0,s2[0]):
        q2New[l,]=n*qF2[l]*f2A[l]**(-1)*fHRU2[l,]
        q2Old[l,]=(1-n)*qF2[l]*f2A[l]**(-1)*fHRU2[l,]+qS2[l]
    
    q2Old=np.nan_to_num(q2Old)
    q2New=np.nan_to_num(q2New)
    
    
    #% P losses (mg/l))
    DRP1=0.0852*pWell-0.3039
    DRP1=np.tile(DRP1,(len(q1New),1))
    
    DRP2=0.0852*pPoor-0.3039
    DRP2=np.tile(DRP2,(len(q2New),1))
                 
    del fHRU1
    del fHRU2
    del f1A
    del f2A
    del qF1
    del qF2
    del qS1
    del qS2
    
    #P load in g
    L1Old=np.sum(q1Old*DRPbase*xres*xres,1)/1000
    L1New=np.sum(q1New*DRP1*xres*xres,1)/1000
    
    del q1Old
    del q1New
    del DRP1
    
    #P load in g    
    L2Old=np.sum(q2Old*DRPbase*xres*xres,1)/1000
    L2New=np.sum(q2New*DRP2*xres*xres,1)/1000
    
    del q2Old
    del q2New
    del DRP2
     
    #Total P load in g
    LTot[:,k]=(L1Old+L2Old+L1New+L2New)
    QTot[:,k]=Q.ravel()
      
    del L1Old
    del L2Old
    del L1New
    del L2New
 
#%% Save results

#Create fast flow map 
h1=h1/nf
h2=h2/nf
hru=np.loadtxt(pathI+'hru.txt',skiprows=6)
m1=np.nonzero(hru==3)
m2=np.nonzero(hru==4)
m3=np.nonzero(hru==2)
m4=np.nonzero(hru==1)
hru[m1]=h1
hru[m2]=h2
hru[m3]=-9999
hru[m4]=-9999

#Write header
header = "ncols     %s\n" % hru.shape[1] 
header += "nrows     %s\n" % hru.shape[0] 
header += "xllcorner   %s\n" % lx
header += "yllcorner    %s\n" % ly
header += "cellsize    %s\n" % xres
header += "NODATA_value    -9999"

#Round results    
resultsPL=pd.DataFrame(LTot.round(3)) #total load (g)
resultsQ=pd.DataFrame(QTot.round(3)) #discharge (mm)

#Save results
resultsQqm=resultsQ*area/1000 # Q in cubicmeter
hourlyC=resultsPL/resultsQqm*1000 #hourly concentrations in mg/m3

np.savetxt(pathO+'hourlyPC.txt',hourlyC,fmt='%1.3f')
np.savetxt(pathO+'hourlyPL.txt',resultsPL,fmt='%1.3f')
np.savetxt(pathO+'hourlyQ.txt',resultsQ,fmt='%1.3f')
np.savetxt(pathO+'hydRisk.asc', hru,header=header,fmt='%1.2f',comments='')

print("Phosphorus module successful!")






