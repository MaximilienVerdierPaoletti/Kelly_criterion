# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 11:41:25 2019

@author: mverdier
"""

try:
    from IPython import get_ipython
    get_ipython().magic('clear')
    get_ipython().magic('reset -f')
except:
    pass

#import h5py #Import necessary module to HDF5 datda extraction
import os #Import module that enable management of paths
import numpy as np
import pandas as pd # enables the use of dataframe
#from scipy import sparse as sp
import matplotlib
import matplotlib.pyplot as plt # Enables plotting of data
import pylab # Enables the use of the 'savefig' function
#from random import randint
import glob #enables to search for specific files in directory and save their names in a list
import re #enables to remove part of strings and look for parts in it
#import sys
import math
import seaborn as sns #to access color palette settings

#Close all open figures
plt.close("all")

#Set working directory
os.chdir('D:/Work/Postdoc/Simulations/FRUITY/Dredge-ups/')

#%%

RS1312C=0.011237
RS1716O=0.00037305
RS1816O=0.002004965
RS2524Mg=0.126597989
RS2624Mg=0.139381904
RS2928Si=0.050775236
RS3028Si=0.033470671



data=pd.read_excel('./Standard 13Cpocket/13CpockstdDU_table.xlsx',header=0)
#data = data.apply(pd.to_numeric, errors='coerce')

Mass=pd.unique(data['Mass'])
Z=pd.unique(data['Metallicity'])


M=input('What mass ?'+str(list(Mass))+'\n')
M=float(M)

dat_Mreduced=data.loc[data['Mass']==M] #Isolate data of the mass of interest
Zrange=pd.unique(dat_Mreduced['Metallicity']) #Extract range of metallicities available for this mass





colors = matplotlib.cm.rainbow(np.linspace(0, 1, len(Zrange))) #Define colors based on metallicity
#dummy_lines = []
#linestyles=['.-','.--']

####################################################  Figures configuration
fig=plt.figure(1)
plt.rcParams['axes.linewidth'] = 1 # Width of borders
plt.rcParams['xtick.major.size'] = 1
plt.rcParams['xtick.major.width'] = 1
plt.rcParams['xtick.minor.size'] = 0.7
plt.rcParams['xtick.minor.width'] = 0.7
manager = plt.get_current_fig_manager() # To open figure in fullscreen
manager.window.showMaximized()


#for b_idx in range(0,2): # Creating dummy lines with no values for additionnal legend
#    dummy_lines.append(plt.plot([],[], c="black", ls = linestyles[b_idx])[0]) 

index=np.array(range(0,data.shape[1])) #Creation of an index to plot the data in the right order


maxL=0
for j in range(0,len(Zrange)):
    print(Zrange[j])
    leg=str(Zrange[j])
    temp_dat=dat_Mreduced.loc[dat_Mreduced['Metallicity']==Zrange[j]]
    temp_dat=temp_dat.dropna(how='all', axis=1)# Remove empty columns
    
    ####################################################  C/O ratio definition
    C12=temp_dat.loc[temp_dat['Isotope']=='C12']
    O16=temp_dat.loc[temp_dat['Isotope']=='O16']
    COrat=C12.iloc[0,5:]/O16.iloc[0,5:]
    
    
    ind=np.max(np.where(COrat.loc[COrat<1]))+6 #Find the TDU event where C/O gets <1 (+5 because of the five first non considered column of C12 and O16 and +1 since index stars at 0)
    
    
    ####################################################  Oxygen ratios 
    O17=temp_dat.loc[temp_dat['Isotope']=='O17']
    O18=temp_dat.loc[temp_dat['Isotope']=='O18']
    O1716=O17.iloc[0,5:ind]/O16.iloc[0,5:ind]
    O1816=O18.iloc[0,5:ind]/O16.iloc[0,5:ind]
    
    
    ####################################################  Mg ratios
    Mg24=temp_dat.loc[temp_dat['Isotope']=='Mg24']
    Mg25=temp_dat.loc[temp_dat['Isotope']=='Mg25']
    Mg26=temp_dat.loc[temp_dat['Isotope']=='Mg26']
    Mg2524=Mg25.iloc[0,5:ind]/Mg24.iloc[0,5:ind]
    Mg2624=Mg26.iloc[0,5:ind]/Mg24.iloc[0,5:ind]
    d25=(Mg2524/RS2524Mg-1)*1000
    d26=(Mg2624/RS2624Mg-1)*1000
    
    
    ####################################################  Si ratios
    Si28=temp_dat.loc[temp_dat['Isotope']=='Si28']
    Si29=temp_dat.loc[temp_dat['Isotope']=='Si29']
    Si30=temp_dat.loc[temp_dat['Isotope']=='Si30']
    Si2928=Si29.iloc[0,5:ind]/Si28.iloc[0,5:ind]
    Si3028=Si30.iloc[0,5:ind]/Si28.iloc[0,5:ind]
    d29=(Si2928/RS2928Si-1)*1000
    d30=(Si3028/RS3028Si-1)*1000  
    

    
    
    
    plt.subplot(221)
    plt.plot(index[5:len(COrat)+5],COrat,'.-',color=colors[j],label=leg)
#    plt.plot(1,'--',color='k')
    plt.xticks(rotation='vertical')
    
    plt.subplot(222)
    plt.plot(index[5:ind],d25,'.-',color=colors[j],label="d25")
    plt.plot(index[5:ind],d26,'.--',color=colors[j],label="d26")
    
    if maxL < len(d25): #Conditional loop to change xticks with the names of maximum TDU events considered
        maxL=len(d25)
        Xlab=list(data.columns[index[5:ind]])
        plt.xticks(range(5,len(Xlab)+5),Xlab,rotation='vertical')
#    plt.xticks(d25,C12.columns[5:-1])
        
    plt.subplot(2,2,3)
    plt.plot(index[5:ind],d29,'.-',color=colors[j])
    plt.plot(index[5:ind],d30,'.--',color=colors[j])
    
    plt.subplot(2,2,4)
    plt.semilogy(index[5:ind],O1716,'.-',color=colors[j])
    plt.semilogy(index[5:ind],O1816,'.--',color=colors[j])
#    plt.plot(index[5:ind],d29/d30,'.-',color=colors[j])
#    plt.plot(index[5:ind],d30,'.--',color=colors[j])
    
    
    
plt.subplot(221)
plt.legend()  
plt.xticks(range(5,len(Xlab)+5),Xlab,rotation='vertical')
plt.subplot(222)
plt.legend()
plt.subplot(223)        
plt.xticks(range(5,len(Xlab)+5),Xlab,rotation='vertical')
plt.subplot(224)
plt.xticks(range(5,len(Xlab)+5),Xlab,rotation='vertical')


#%%
'''%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Plot basic definition   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'''


#Extracting measurements
meas=pd.read_excel('D:/Work/Postdoc/NanoSIMS/Presolar silicates/Acfer094/Table.xlsx',header=0)


fig=plt.figure(2)


plt.rcParams['axes.linewidth'] = 1 # Width of borders
plt.rcParams['xtick.major.size'] = 1
plt.rcParams['xtick.major.width'] = 1
plt.rcParams['xtick.minor.size'] = 0.7
plt.rcParams['xtick.minor.width'] = 0.7
manager = plt.get_current_fig_manager() # To open figure in fullscreen
manager.window.showMaximized()




############################################ Figure 1
FT=15
loc_iso='δ25Mg (‰)'


plt.subplot(221)
plt.plot(meas.loc[(meas[str(loc_iso)]<0) & (meas['δ26Mg (‰)']>0),'δ26Mg (‰)'],meas.loc[(meas[str(loc_iso)]<0) & (meas['δ26Mg (‰)']>0),'δ25Mg (‰)'],marker='o', markerfacecolor='r',linestyle='None',zorder=10,label='_nolegend_')
plt.errorbar(meas['δ26Mg (‰)'],meas['δ25Mg (‰)'],yerr=meas['σ.3']*2,xerr=meas['σ.2']*2,fmt='o',color='grey',label='measurements',zorder=5)
plt.xlim([-200,1100])
plt.ylim([-200,2500])
plt.vlines(0,round(plt.gca().get_ylim()[0])-100,math.ceil(plt.gca().get_ylim()[1]),'k',linestyle='dashed')
plt.hlines(0,round(plt.gca().get_xlim()[0])-100,math.ceil(plt.gca().get_xlim()[1]),'k',linestyle='dashed')
plt.xlabel(u'$\delta^{26}$Mg (‰)',fontsize=FT)
plt.ylabel(u'$\delta^{25}$Mg (‰)',fontsize=FT)

plt.subplot(222)
plt.xscale('log', nonposx='clip')
plt.yscale('log', nonposy='clip')
plt.plot(meas.loc[(meas[str(loc_iso)]<0) & (meas['δ26Mg (‰)']>0),'18O/16O(x10-4)']*(1e-4),meas.loc[(meas[str(loc_iso)]<0) & (meas['δ26Mg (‰)']>0),'17O/16O(x10-4)']*(1e-4),marker='o', markerfacecolor='r',linestyle='None',zorder=10,label='_nolegend_')
plt.errorbar(meas['18O/16O(x10-4)']*(1e-4),meas['17O/16O(x10-4)']*(1e-4),yerr=meas['σ']*(1e-4),xerr=meas['σ.1']*(1e-4),fmt='o',color='grey',label='measurements',zorder=5)
plt.vlines(RS1816O,round(plt.gca().get_ylim()[0]),math.ceil(plt.gca().get_ylim()[1]),'k',linestyle='dashed')
plt.hlines(RS1716O,round(plt.gca().get_xlim()[0]),math.ceil(plt.gca().get_xlim()[1]),'k',linestyle='dashed')
plt.xlim([1e-3,3e-3])
plt.ylim([1e-4,1e-2])
plt.xlabel(u'$^{18}$O/$^{16}$O',fontsize=FT)
plt.ylabel(u'$^{17}$O/$^{16}$O',fontsize=FT)


plt.subplot(223,aspect='equal')
plt.xscale('log', nonposx='clip')
plt.plot(meas.loc[(meas[str(loc_iso)]<0) & (meas['δ26Mg (‰)']>0),'18O/16O(x10-4)']*(1e-4),meas.loc[(meas[str(loc_iso)]<0) & (meas['δ26Mg (‰)']>0),'δ25Mg (‰)'],marker='o', markerfacecolor='r',linestyle='None',zorder=10,label='_nolegend_')
plt.errorbar(meas['18O/16O(x10-4)']*(1e-4),meas['δ25Mg (‰)'],yerr=meas['σ.3'],xerr=meas['σ.1']*(1e-4),fmt='o',color='grey',label='measurements',zorder=5)
plt.xlim([1e-3,3e-3])
plt.ylim([-200,2500])
plt.vlines(RS1816O,-200,math.ceil(plt.gca().get_ylim()[1]),'k',linestyle='dashed')
plt.hlines(0,round(plt.gca().get_xlim()[0]),math.ceil(plt.gca().get_xlim()[1]),'k',linestyle='dashed')
plt.xlabel(u'$^{18}$O/$^{16}$O',fontsize=FT)
plt.ylabel(u'$\delta^{25}$Mg (‰)',fontsize=FT)



plt.subplot(224)
plt.plot(meas.loc[(meas[str(loc_iso)]<0) & (meas['δ26Mg (‰)']>0),'δ30Si (‰)'],meas.loc[(meas[str(loc_iso)]<0) & (meas['δ26Mg (‰)']>0),'δ29Si (‰)'],marker='o', markerfacecolor='r',linestyle='None',zorder=10,label='_nolegend_')
plt.errorbar(meas['δ30Si (‰)'],meas['δ29Si (‰)'],yerr=meas['σ.4'],xerr=meas['σ.5'],fmt='o',color='grey',label='measurements',zorder=5)
plt.vlines(0,round(plt.gca().get_ylim()[0])-100,math.ceil(plt.gca().get_ylim()[1])+100,'k',linestyle='dashed')
plt.hlines(0,round(plt.gca().get_xlim()[0])-100,math.ceil(plt.gca().get_xlim()[1])+100,'k',linestyle='dashed')
plt.xlim([-50,200])
plt.ylim([-20,250])
plt.xlabel(u'$\delta^{30}$Si (‰)',fontsize=FT)
plt.ylabel(u'$\delta^{29}$Si (‰)',fontsize=FT)



#%%
'''%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Plot basic definition   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'''

M=2
Z=0.003
TDU_num='TDU_1'

mod=data.loc[data['Mass']==M] #Isolate data of the mass of interest
mod=mod.loc[mod['Metallicity']==Z]


dilu=np.arange(0,1.05,0.05) #Percentage of dilution array from 0 to 100%
R17d=np.zeros(len(dilu)-1) #Allocate space to d25 array
R18d=np.zeros(len(dilu)-1) #Allocate space to d26 array
d25d=np.zeros(len(dilu)-1) #Allocate space to d25 array
d26d=np.zeros(len(dilu)-1) #Allocate space to d26 array
d29d=np.zeros(len(dilu)-1) #Allocate space to d25 array
d30d=np.zeros(len(dilu)-1) #Allocate space to d26 array

O16=float(mod.loc[mod['Isotope']=='O16'][TDU_num])
O17=float(mod.loc[mod['Isotope']=='O17'][TDU_num])
O18=float(mod.loc[mod['Isotope']=='O18'][TDU_num])
Mg24=float(mod.loc[mod['Isotope']=='Mg24'][TDU_num])
Mg25=float(mod.loc[mod['Isotope']=='Mg25'][TDU_num])
Mg26=float(mod.loc[mod['Isotope']=='Mg26'][TDU_num])
Si28=float(mod.loc[mod['Isotope']=='Si28'][TDU_num])
Si29=float(mod.loc[mod['Isotope']=='Si29'][TDU_num])
Si30=float(mod.loc[mod['Isotope']=='Si30'][TDU_num])

for i in range(0,len(dilu)-1):
    R17d[i]=(O17/O16)*dilu[i]+(1-dilu[i])*RS1716O
    R18d[i]=(O18/O16)*dilu[i]+(1-dilu[i])*RS1816O
    
    Mg2524=Mg25/Mg24
    Mg2624=Mg26/Mg24
    d25d[i]=((Mg2524*dilu[i])/RS2524Mg-dilu[i])*1000
    d26d[i]=((Mg2624*dilu[i])/RS2624Mg-dilu[i])*1000
    
    Si2928=Si29/Si28
    Si3028=Si30/Si28
    d29d[i]=((Si2928*dilu[i])/RS2928Si-dilu[i])*1000
    d30d[i]=((Si3028*dilu[i])/RS3028Si-dilu[i])*1000 
    
plt.subplot(221)
plt.plot(d26d,d25d,'.-b')

plt.subplot(222)
plt.loglog(R18d,R17d,'.-b')

plt.subplot(223,aspect='equal')
plt.semilogx(R18d,d25d,'.-b')

plt.subplot(224)
plt.plot(d30d,d29d,'.-b')




