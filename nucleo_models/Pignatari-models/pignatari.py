# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 16:12:54 2018

@author: mverdier
"""


from IPython import get_ipython
if '__IPYTHON__' in globals():
    get_ipython().magic('reset -sf') #Give access to the command to clear all the variables at once


import h5py #Import necessary module to HDF5 datda extraction
import os #Import module that enable management of paths
import numpy as np
import pandas as pd # enables the use of dataframe
from scipy import sparse as sp
import matplotlib.pyplot as plt # Enables plotting of data
import pylab # Enables the use of the 'savefig' function
from random import randint
import glob #enables to search for specific files in directory and save their names in a list
import re #enables to remove part of strings


os.chdir('D:/Work/Programmation/Presolar grains/Pignatari-models/')


plt.close("all")


#%%
'''%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Variable declarations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'''

FT=15
species=['C_', 'Mg_', 'Al_', 'Si_']
CSizonecoord=[6.83,7.04]
models=['m25_exp_d','m25_exp_d_all_H10','m25_exp_d_all_H500','m25_exp_T_all','m25_exp_T_all_H10','m25_exp_T_all_H500']

print("\n"*100)
m=input('Which mass-loss model (type numer starting from 0): '+str(models)+'\n')
m=int(m)



#%%
'''%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Extracting data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'''


data=pd.read_table(str(models[m])+'.txt',sep=" ", header=0) #Reads the txt considering space as delimiter and the first row (starts at 0) as the names of columns
pd.DataFrame.describe(data) #To get a summary of the data


colors = []

for i in range(len(list(data))-1): #Loops over the number of species in the data
    colors.append('%06X' % randint(0, 0xFFFFFF))
    
#colnames=[]
#for i in range(0,len(species)):
#    colnames.append([col for col in data.columns if names[i] in col])
    

#%% 
    '''%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Display of species of interest  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'''
    
    
fig=plt.figure()
fig.suptitle('Supernovae model with H ingestion '+str(models[m])+' Pignatari et al. (2015)', fontsize=16)
#plot basic definition
plt.rcParams['axes.linewidth'] = 1 # Width of borders
plt.rcParams['xtick.major.size'] = 1
plt.rcParams['xtick.major.width'] = 1
plt.rcParams['xtick.minor.size'] = 0.7
plt.rcParams['xtick.minor.width'] = 0.7
manager = plt.get_current_fig_manager() # To open figure in fullscreen
manager.window.showMaximized()


for j in range(0,len(species)):
    plt.subplot(2,2,j+1)
    a=np.shape(data.filter(regex=species[j])); #Isolate data of species of interest
    for i in range(0,a[1]):# Loop over the isolated data
        plt.semilogy(data.mass,data.filter(regex=species[j]).iloc[:,i]) 
    plt.legend()
    plt.axvspan(CSizonecoord[0],CSizonecoord[1],facecolor='grey',alpha=0.5) #color in grey the region of the C/Si region
    plt.text(6.935,-0.5,"C/Si",fontsize=12,horizontalalignment='center')
#    plt.text(7.32,-0.5,"He/C",fontsize=12,horizontalalignment='center')
    plt.xlabel('Mass coordinates ($M_{\odot}$)',fontsize=FT)
    plt.ylabel('Mass fraction',fontsize=FT)
    


#%%
#
#C12_13=data.C_12/data.C_13
#N14_15=data.N_14/data.N_15
#
##ind_HeC=np.where(data.mass>7.04)
#plt.figure()
##plt.loglog(C12_13[ind_HeC[0]], N14_15[ind_HeC[0]])
#plt.loglog(C12_13, N14_15)
#plt.axis([ 10**(-3), 10**4,1, 10**5])
#
#plt.figure()
#plt.plot(data.mass,np.log10(C12_13))
#plt.axvspan(6.83,7.04,facecolor='grey',alpha=0.5) #color in grey the region of the C/Si region
#plt.text(6.935,10,"C/Si",fontsize=12,horizontalalignment='center')
#plt.text(7.32,10,"He/C",fontsize=12,horizontalalignment='center')

#%% Ti isotopes in SN

#files_names=glob.glob("m*.txt")
#separator=[" ","\t"," ", " ","\t","\t"]
#
#f, ax = plt.subplots(2, 3)
#
#for i in range(1,3):
#    print(i)
#    for j in range(1,4):
#        print(j)
#
##           data=pd.read_table(files_names[i],sep=" ", header=0)
##           ind_HeC=np.where((data.mass>6.81) & (data.mass<9.23))
#    
#    
##    C12_13=data.C_12[ind_HeC[0]]/data.C_13[ind_HeC[0]]
##    N14_15=data.N_14[ind_HeC[0]]/data.N_15[ind_HeC[0]]
#    
#            
#            
#        if i==1:
##           lab=lb
#            lb=re.sub('.txt', '', files_names[j-i])
#            data=pd.read_table(files_names[j-i],sep=separator[j-i], header=0)
#        else:
##           lab=[lab,lb]
#            lb=re.sub('.txt', '', files_names[j+i])
#            data=pd.read_table(files_names[j+i],sep=separator[j+i], header=0)
#                
#        lb=re.sub('_exp', '', lb)
##        ind_HeC=np.where(data.mass>6.82) #From the C/Si zone through the O/nova zone until the end of the He/C zone
#        ax[i-1,j-1].loglog(data.mass,data.Ti_46)
#        ax[i-1,j-1].loglog(data.mass,data.Ti_47)
#        ax[i-1,j-1].loglog(data.mass,data.Ti_49)
#        ax[i-1,j-1].loglog(data.mass,data.Ti_50)
#        ax[i-1,j-1].axvspan(6.83,7.04,facecolor='grey',alpha=0.5)
#        ax[i-1,j-1].set_title(lb)
#        ax[i-1,j-1].legend(loc='lower left')
#            
#    
##plt.figlegend(lines, loc = 'lower center', ncol=5, labelspacing=0. )
##plt.plot([1E-6,1E6],[440.91,440.91],'k--',lw=2) #Marty2011 Solar 14N/15N
##plt.plot([89.9,89.9],[1E-6,1E6],'k--',lw=2) #Solar 12C/13C
##plt.axis([ 10**(-3), 10**4,1, 10**5])



