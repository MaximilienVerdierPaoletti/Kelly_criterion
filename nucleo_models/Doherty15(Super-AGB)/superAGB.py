# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 12:15:54 2019

@author: mverdier
"""

#Super-AGB stars models based on Doherty2015 paper II

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
os.chdir('D:/Work/Programmation/Presolar grains/')


#%%

############# Figure variables
linecolors='w'
fontcolor=linecolors
FT=30
MS_ini=10
MS=12
generalcolor=linecolors

plt.rcParams['axes.linewidth'] = 2 # Width of borders
#plt.rcParams['axes.facecolor'] = 'k'
plt.rcParams['axes.edgecolor'] = generalcolor
plt.rcParams['lines.color'] = generalcolor
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['xtick.major.size'] = 8
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['xtick.minor.size'] = 4
plt.rcParams['xtick.minor.width'] = 2
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['xtick.color'] = generalcolor
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['ytick.major.size'] = 8
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['ytick.minor.size'] = 4
plt.rcParams['ytick.minor.width'] = 2
plt.rcParams['ytick.labelsize'] = 20
plt.rcParams['ytick.color'] = generalcolor


#%%
'''%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Extraction of data   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'''




#Extracting measurements
data=pd.read_excel('D:/Work/NanoSIMS Data/Presolar grains/Presolar Silicates/Acfer094/Table.xlsx',header=0)





print("\n"*100)
A=input('AGB-star, Super-AGB-stars or both ?(0,1 or 2)\n')
A=int(A)
if A == 1:
    models_names=['VW93','VW-M','testcases']

    
    #Extract models data
    pd.set_option("display.precision", 8)
    Mini=pd.read_excel('./Doherty15(Super-AGB)/InitialCompositions.xlsx',header=0)
    
    
    m=input('Which mass-loss model: '+str(models_names)+'(0 - 1 - 2)\n')
    m=int(m)
    mod=pd.read_excel('./Doherty15(Super-AGB)/reformatted'+str(models_names[m])+'.xlsx',header=0) #"E-" columns correspond to yields corresponding to extrapolated Thermal Pulses (TP) after convergence of the main model
    
    
    starmass=mod.Mass.unique()
    starmet=mod.Metallicity.unique()
    
    
    
elif A==0: #Fruity models
    models_names=['Fruity']
    m=0
    mod=pd.read_excel('./FRUITY/fruity_totalyields.xlsx',header=0)
    starmass=mod.Mass.unique()
    starmet=mod.Metallicity.unique()
    
    
else:
    print("Unavailable at the moment")
            

#%%
'''%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Settings of constant and memory allocation   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'''

#Diverse
loc_iso='δ25Mg (‰)'


#Solar ratios (Lodders 2003)
RS1312C=0.011237
RS1716O=0.00037305
RS1816O=0.002004965
RS2524Mg=0.126597989
RS2624Mg=0.139381904
RS2928Si=0.050775236
RS3028Si=0.033470671


dilu=np.arange(0,1.05,0.05) #Percentage of dilution array from 0 to 100%
#Memory allocation
R1716O=np.zeros(len(dilu)-1)
R1816O=np.zeros(len(dilu)-1)
ER1716O=np.zeros(len(dilu)-1)
ER1816O=np.zeros(len(dilu)-1)
        
        
d25=np.zeros(len(dilu)-1) #Allocate space to d25 array
d26=np.zeros(len(dilu)-1) #Allocate space to d26 array
Ed25=np.zeros(len(dilu)-1) #Allocate space to d25 array
Ed26=np.zeros(len(dilu)-1) #Allocate space to d26 array
        
d29=np.zeros(len(dilu)-1) #Allocate space to d25 array
d30=np.zeros(len(dilu)-1) #Allocate space to d26 array
Ed29=np.zeros(len(dilu)-1) #Allocate space to d25 array
Ed30=np.zeros(len(dilu)-1) #Allocate space to d26 array



#Ask the user for the desired dataset
i=input('Vary star masses or metallicities? (M or Z) \n')
Var=str(i)



#%%
'''%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Function definitions   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'''




def funcgroup(O18,O17):
   if len(O18) != len(O17):
       print("Error, argument must be of same length")
       return
   else:
#       group=np.zeros(len(O18)-1)
       group=list()
       for i in range(0,len(O18)):
           if (O18[i] <0.0008 and O17[i]>RS1716O):
#               group[i]=int(2)
               group.append(2)
           elif (O17[i]<RS1716O and O18[i]<0.00198):
#               group[i]=int(3)
               group.append(3)
           elif O18[i]>0.003125:
#               group[i]=int(4)
               group.append(4)
           else:
#               group[i]=int(1)
               group.append(1)
       return group



#%%
'''%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Group definition   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'''


gr=funcgroup(data['18O/16O(x10-4)']*(1e-4),data['17O/16O(x10-4)']*(1e-4))
#data['Group']=gr

    
#%%
'''%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Plot basic definition   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'''


FT=15
dummy_lines = []
linestyles=['-','--']




fig, [[ax1,ax2],[ax3,ax4]]=plt.subplots(2,2,figsize=(10,10))

#colorscale
if Var == 'M':
    colors = matplotlib.cm.rainbow(np.linspace(0, 1, len(starmass)))
    L=len(starmass)-1 # We do not consider mass 6.5
    Z=input('What metallicty ?'+str(list(np.unique(starmet)))+'\n')
    Z=float(Z)
    fig.suptitle('Star mass dependance, mass-loss model '+str(models_names[m])+', Z='+str(Z), fontsize=16)
    
    if A==1:
        #Initial Ratios
    #    RI1312C=Mini.loc[Mini['Species']=='c13','Z='+str(Z)].item()/Mini.loc[Mini['Species']=='c12','Z='+str(Z)].item()
        RI1716O=Mini.loc[Mini['Species']=='o17','Z='+str(Z)].item()/Mini.loc[Mini['Species']=='o16','Z='+str(Z)].item()
        RI1816O=Mini.loc[Mini['Species']=='o18','Z='+str(Z)].item()/Mini.loc[Mini['Species']=='o16','Z='+str(Z)].item()
        RI2524Mg=Mini.loc[Mini['Species']=='mg25','Z='+str(Z)].item()/Mini.loc[Mini['Species']=='mg24','Z='+str(Z)].item()
        RI2624Mg=Mini.loc[Mini['Species']=='mg26','Z='+str(Z)].item()/Mini.loc[Mini['Species']=='mg24','Z='+str(Z)].item()
        RI2928Si=Mini.loc[Mini['Species']=='si29','Z='+str(Z)].item()/Mini.loc[Mini['Species']=='si28','Z='+str(Z)].item()
        RI3028Si=Mini.loc[Mini['Species']=='si30','Z='+str(Z)].item()/Mini.loc[Mini['Species']=='si28','Z='+str(Z)].item() 
else:
    colors = matplotlib.cm.rainbow(np.linspace(0, 1, len(starmet)))
    L=len(starmet)
    mass=input('What mass ?'+str(list(np.unique(starmass)))+'\n')
    mass=float(mass)
    fig.suptitle('Star metallicity dependance, mass-loss model '+str(models_names[m])+ ', M='+str(mass)+'$M_{\odot}$', fontsize=16)




#plot basic definition
# plt.rcParams['axes.linewidth'] = 1 # Width of borders
# plt.rcParams['xtick.major.size'] = 1
# plt.rcParams['xtick.major.width'] = 1
# plt.rcParams['xtick.minor.size'] = 0.7
# plt.rcParams['xtick.minor.width'] = 0.7
manager = plt.get_current_fig_manager() # To open figure in fullscreen
manager.window.showMaximized()


ax1.plot(data.loc[data[str(loc_iso)]>1900,'δ26Mg (‰)'],data.loc[data[str(loc_iso)]>1900,'δ25Mg (‰)'],marker='o', markerfacecolor='r',linestyle='None',zorder=10,label='_nolegend_')
ax1.errorbar(data['δ26Mg (‰)'],data['δ25Mg (‰)'],yerr=data['σ.3']*2,xerr=data['σ.2']*2,fmt='o',color='cornflowerblue',label='data',zorder=5)
ax1.set_xlim([-200,1100])
ax1.set_ylim([-200,2500])
ax1.vlines(0,round(ax1.get_ylim()[0])-100,math.ceil(ax1.get_ylim()[1]),'k',linestyle='dashed')
ax1.hlines(0,round(ax1.get_xlim()[0])-100,math.ceil(ax1.get_xlim()[1]),'k',linestyle='dashed')
ax1.set_xlabel(u'$\delta^{26}$Mg (‰)',fontsize=FT)
ax1.set_ylabel(u'$\delta^{25}$Mg (‰)',fontsize=FT)

ax2.set_xscale('log', nonposx='clip')
ax2.set_yscale('log', nonposy='clip')
ax2.plot(data.loc[data[str(loc_iso)]>1900,'18O/16O(x10-4)']*(1e-4),data.loc[data[str(loc_iso)]>1900,'17O/16O(x10-4)']*(1e-4),marker='o', markerfacecolor='r',linestyle='None',zorder=10,label='_nolegend_')
ax2.errorbar(data['18O/16O(x10-4)']*(1e-4),data['17O/16O(x10-4)']*(1e-4),yerr=data['σ']*(1e-4),xerr=data['σ.1']*(1e-4),fmt='o',color='cornflowerblue',label='data',zorder=5)
ax2.vlines(RS1816O,round(ax2.get_ylim()[0]),math.ceil(ax2.get_ylim()[1]),'k',linestyle='dashed')
ax2.hlines(RS1716O,round(ax2.get_xlim()[0]),math.ceil(ax2.get_xlim()[1]),'k',linestyle='dashed')
ax2.set_xlim([1e-3,3e-3])
ax2.set_ylim([1e-4,1e-2])
ax2.set_xlabel(u'$^{18}$O/$^{16}$O',fontsize=FT)
ax2.set_ylabel(u'$^{17}$O/$^{16}$O',fontsize=FT)


ax3.set_xscale('log', nonposx='clip')
ax3.plot(data.loc[data[str(loc_iso)]>1900,'18O/16O(x10-4)']*(1e-4),data.loc[data[str(loc_iso)]>1900,'δ25Mg (‰)'],marker='o', markerfacecolor='r',linestyle='None',zorder=10,label='_nolegend_')
ax3.errorbar(data['18O/16O(x10-4)']*(1e-4),data['δ25Mg (‰)'],yerr=data['σ.3'],xerr=data['σ.1']*(1e-4),fmt='o',color='cornflowerblue',label='data',zorder=5)
ax3.set_xlim([1e-3,3e-3])
ax3.set_ylim([-200,2500])
ax3.vlines(RS1816O,-200,math.ceil(ax3.get_ylim()[1]),'k',linestyle='dashed')
ax3.hlines(0,round(ax3.get_xlim()[0]),math.ceil(ax3.get_xlim()[1]),'k',linestyle='dashed')
ax3.set_xlabel(u'$^{18}$O/$^{16}$O',fontsize=FT)
ax3.set_ylabel(u'$\delta^{25}$Mg (‰)',fontsize=FT)



ax4.plot(data.loc[data[str(loc_iso)]>1900,'δ30Si (‰)'],data.loc[data[str(loc_iso)]>1900,'δ29Si (‰)'],marker='o', markerfacecolor='r',linestyle='None',zorder=10,label='_nolegend_')
ax4.errorbar(data['δ30Si (‰)'],data['δ29Si (‰)'],yerr=data['σ.4'],xerr=data['σ.5'],fmt='o',color='cornflowerblue',label='data',zorder=5)
ax4.vlines(0,round(ax4.get_ylim()[0])-100,math.ceil(ax4.get_ylim()[1])+100,'k',linestyle='dashed')
ax4.hlines(0,round(ax4.get_xlim()[0])-100,math.ceil(ax4.get_xlim()[1])+100,'k',linestyle='dashed')
ax4.set_xlim([-50,200])
ax4.set_ylim([-20,250])
ax4.set_xlabel(u'$\delta^{30}$Si (‰)',fontsize=FT)
ax4.set_ylabel(u'$\delta^{29}$Si (‰)',fontsize=FT)





fig_evol=plt.figure(2)
#%%


fig2= plt.figure(figsize=(10,10))

plt.plot(data.loc[data[str(loc_iso)]>1900,'δ26Mg (‰)'],data.loc[data[str(loc_iso)]>1900,'δ25Mg (‰)'],marker='o', markerfacecolor='r',linestyle='None',markersize=MS,zorder=10,label='_nolegend_')
plt.errorbar(data['δ26Mg (‰)'],data['δ25Mg (‰)'],yerr=data['σ.3']*2,xerr=data['σ.2']*2,fmt='o',color='cornflowerblue',markersize=MS,label='data',zorder=5)
plt.axis([-200,1100,-200,2500])
plt.vlines(0,round(plt.gca().get_ylim()[0])-100,math.ceil(plt.gca().get_ylim()[1]),'w',linestyle='dashed')
plt.hlines(0,round(plt.gca().get_xlim()[0])-100,math.ceil(plt.gca().get_xlim()[1]),'w',linestyle='dashed')
plt.xlabel(u'$\delta^{26}$Mg (‰)',fontsize=FT,color='w')
plt.ylabel(u'$\delta^{25}$Mg (‰)',fontsize=FT,color='w')


# fig2.savefig('G:/Mon Drive/Mg.png', transparent=True)






#%%
'''%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Loops %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'''

for j in range(0,L):

    if Var== 'M':    
        mass=starmass[j]
        leg=str(starmass[j])+'$M_{\odot}$'
        Xaxis=mass
    else:
        Z=starmet[j]
        leg='Z='+str(starmet[j])
        Xaxis=Z
        if A==1:    
            #Initial ratios
    #        RI1312C=Mini.loc[Mini['Species']=='c13','Z='+str(Z)].item()/Mini.loc[Mini['Species']=='c12','Z='+str(Z)].item()
            RI1716O=Mini.loc[Mini['Species']=='o17','Z='+str(Z)].item()/Mini.loc[Mini['Species']=='o16','Z='+str(Z)].item()
            RI1816O=Mini.loc[Mini['Species']=='o18','Z='+str(Z)].item()/Mini.loc[Mini['Species']=='o16','Z='+str(Z)].item()
            RI2524Mg=Mini.loc[Mini['Species']=='mg25','Z='+str(Z)].item()/Mini.loc[Mini['Species']=='mg24','Z='+str(Z)].item()
            RI2624Mg=Mini.loc[Mini['Species']=='mg26','Z='+str(Z)].item()/Mini.loc[Mini['Species']=='mg24','Z='+str(Z)].item()
            RI2928Si=Mini.loc[Mini['Species']=='si29','Z='+str(Z)].item()/Mini.loc[Mini['Species']=='si28','Z='+str(Z)].item()
            RI3028Si=Mini.loc[Mini['Species']=='si30','Z='+str(Z)].item()/Mini.loc[Mini['Species']=='si28','Z='+str(Z)].item()           
        
    
    #Locate data of interest
    mod_reduced=mod.loc[(mod['Mass'] == mass) & (mod['Metallicity'] == Z)]
    if mod_reduced.empty==True: # If no match is found (i.e. mod_reduced is empty) then skip to next iteration
        continue
    
    
    if A==1: #Super-AGB models
        for b_idx in range(0,2): # Creating dummy lines with no values for additionnal legend
            dummy_lines.append(plt.plot([],[], c="black", ls = linestyles[b_idx])[0]) 
    
        
        #Extract Prod factors (PF) and Extrapolated prod factors
        Y16=mod_reduced.loc[mod_reduced['Species']=='o16','ProdFact'].item()
        Y17=mod_reduced.loc[mod_reduced['Species']=='o17','ProdFact'].item()
        Y18=mod_reduced.loc[mod_reduced['Species']=='o18','ProdFact'].item()
        EY16=mod_reduced.loc[mod_reduced['Species']=='o16','E-ProdFact'].item()
        EY17=mod_reduced.loc[mod_reduced['Species']=='o17','E-ProdFact'].item()
        EY18=mod_reduced.loc[mod_reduced['Species']=='o18','E-ProdFact'].item()
        
        
        Y24=mod_reduced.loc[mod_reduced['Species']=='mg24','ProdFact'].item()
        Y25=mod_reduced.loc[mod_reduced['Species']=='mg25','ProdFact'].item()
        Y26=mod_reduced.loc[mod_reduced['Species']=='mg26','ProdFact'].item()
        EY24=mod_reduced.loc[mod_reduced['Species']=='mg24','E-ProdFact'].item()
        EY25=mod_reduced.loc[mod_reduced['Species']=='mg25','E-ProdFact'].item()
        EY26=mod_reduced.loc[mod_reduced['Species']=='mg26','E-ProdFact'].item()
        
        Y28=mod_reduced.loc[mod_reduced['Species']=='si28','ProdFact'].item()
        Y29=mod_reduced.loc[mod_reduced['Species']=='si29','ProdFact'].item()
        Y30=mod_reduced.loc[mod_reduced['Species']=='si30','ProdFact'].item()
        EY28=mod_reduced.loc[mod_reduced['Species']=='si28','E-ProdFact'].item()
        EY29=mod_reduced.loc[mod_reduced['Species']=='si29','E-ProdFact'].item()
        EY30=mod_reduced.loc[mod_reduced['Species']=='si30','E-ProdFact'].item()
        
    
    
        #loop calculating diluted delta values
        for i in range(0,len(dilu)-1):
            #delta=((solar*(10^PFnumerator/10^PFdenominator)*dilution+(1-dilution)*solar)/solar-1)*1000
            R1716O[i]=RI1716O*(10**Y17)/(10**Y16)*dilu[i]+(1-dilu[i])*RS1716O
            R1816O[i]=RI1816O*(10**Y18)/(10**Y16)*dilu[i]+(1-dilu[i])*RS1816O
            ER1716O[i]=RI1716O*(10**EY17)/(10**EY16)*dilu[i]+(1-dilu[i])*RS1716O
            ER1816O[i]=RI1816O*(10**EY18)/(10**EY16)*dilu[i]+(1-dilu[i])*RS1816O
            
            d25[i]=(((RI2524Mg*(10**Y25)/(10**Y24)*dilu[i])/RS2524Mg)-dilu[i])*1000
            d26[i]=(((RI2624Mg*(10**Y26)/(10**Y24)*dilu[i])/RS2624Mg)-dilu[i])*1000
            Ed25[i]=(((RI2524Mg*(10**EY25)/(10**EY24)*dilu[i])/RS2524Mg)-dilu[i])*1000
            Ed26[i]=(((RI2624Mg*(10**EY26)/(10**EY24)*dilu[i])/RS2624Mg)-dilu[i])*1000
            
            d29[i]=(((RI2928Si*(10**Y29)/(10**Y28)*dilu[i])/RS2928Si)-dilu[i])*1000
            d30[i]=(((RI3028Si*(10**Y30)/(10**Y28)*dilu[i])/RS3028Si)-dilu[i])*1000
            Ed29[i]=(((RI2928Si*(10**EY29)/(10**EY28)*dilu[i])/RS2928Si)-dilu[i])*1000
            Ed30[i]=(((RI3028Si*(10**EY30)/(10**EY28)*dilu[i])/RS3028Si)-dilu[i])*1000
            
    elif A==0:
        for i in range(0,len(dilu)-1):
            R1716O[i]=(mod_reduced.o17/mod_reduced.o16)*dilu[i]+(1-dilu[i])*RS1716O
            R1816O[i]=(mod_reduced.o18/mod_reduced.o16)*dilu[i]+(1-dilu[i])*RS1816O
          
            
            d25[i]=((mod_reduced.mg25/mod_reduced.mg24*dilu[i]/RS2524Mg)-dilu[i])*1000
            d26[i]=((mod_reduced.mg26/mod_reduced.mg24*dilu[i]/RS2624Mg)-dilu[i])*1000
            
            
            d29[i]=((mod_reduced.si29/mod_reduced.si28*dilu[i]/RS2928Si)-dilu[i])*1000
            d30[i]=((mod_reduced.si30/mod_reduced.si28*dilu[i]/RS3028Si)-dilu[i])*1000
            
        
    
       
        
        
              
        
    #plots
    plt.figure(1)
    ax1.plot(d26,d25,'.-',color=colors[j],label=leg)


    ax2.loglog(R1816O,R1716O,'.-',color=colors[j])

    

    ax3.semilogx(R1816O,d25,'.-',color=colors[j])

    

    ax4.plot(d30,d29,'.-',color=colors[j])

    
    if A==1:      
        ax1.plot(Ed26,Ed25,'--',color=colors[j])
        
        ax2.loglog(ER1816O,ER1716O,'--',color=colors[j])

        ax3.semilogx(ER1816O,Ed25,'--',color=colors[j])

        ax4.plot(Ed30,Ed29,'--',color=colors[j])
        
        
    plt.figure(2)
    plt.subplot(221)
    plt.plot(Xaxis,d25[0],'.-',color='r',label=leg)
    plt.plot(Xaxis,d26[0],'.-',color='b',label=leg)

    plt.subplot(222)
    plt.loglog(Xaxis,R1716O[0],'.-',color='r')
    plt.loglog(Xaxis,R1816O[0],'.-',color='b')
    
    plt.subplot(223)
    plt.semilogx(Xaxis,d29[0],'.-',color='r')
    plt.semilogx(Xaxis,d30[0],'.-',color='b')
    
    
    
    
  
plt.figure(1)
legend1=ax1.legend(loc=2)

if A==1:
    legend2=ax1.legend()
    legend2=ax1.legend([dummy_lines[i] for i in [0,1]], ["Converging models", "Converging models + Extrapolated pulses"], loc=1)
    ax1.add_artist(legend1)


    


#%%


M=7
Z=0.008

solar_abundances=pd.read_excel('D:/Work/Abondances.xlsx',sheet_name='isotope ratios',header=0)



# AGB models
models_names=['Fruity']
m=0
AGB=pd.read_excel('./FRUITY/fruity_totalyields.xlsx',header=0)
AGB_mass=AGB.Mass.unique()
AGB_met=AGB.Metallicity.unique()

# S-AGB models. Works only if S-AGB models were selected in the beginning
SAGB=mod

SAGB_mass=SAGB.Mass.unique()
SAGB_met=SAGB.Metallicity.unique()



# Figure definition
fig3, axs= plt.subplots(2,2,figsize=(10,10))
axs=axs.ravel()

met_common=np.sort(list(set(AGB_met) & set(SAGB.Metallicity)))

for h in range(0,len(met_common)):
    Z=met_common[h]
    SAGB_reduced=SAGB.loc[SAGB.Metallicity==Z]
    SAGB_species=SAGB_reduced.Species.unique()
    Y_species=[SAGB_reduced.loc[SAGB_reduced.Species==i,'ProdFact'] for i in SAGB_species]
    
    mass_iso=[int(''.join(list(filter(str.isdigit, s)))) for s in SAGB_species if len(s)>1] #extract mass of elements available in Super AGB models
    letters=[str(''.join(list(filter(str.isalpha, s)))) for s in SAGB_species if len(s)>1]
    SAGBelem_mass=[str(str(mass_iso[i])+letters[i]) for i in range(0,len(letters))]
    
    iso_ratio=solar_abundances[solar_abundances.stack().str.contains('|'.join(SAGBelem_mass),case=False).any(level=0)]
    x=np.linspace(1,len(iso_ratio)+1)
    colors = matplotlib.cm.rainbow(np.linspace(0, 1, len(x)))
    encountered_ratio=[]
    tickposition=[]
    
    for j in SAGB_mass:
        SAGB_rereduced=SAGB_reduced.loc[SAGB_reduced.Mass==j]
    
        for i in range(0,len(iso_ratio)):
            split_string = iso_ratio.ratio_name.iloc[i].split("/", 1)
            major=split_string[1]
            minor=split_string[0]
            ratio_solar=iso_ratio.ratio.iloc[i].item()
            
            major_search=str(str(''.join(list(filter(str.isalpha, major))))+str(''.join(list(filter(str.isdigit, major)))))
            minor_search=str(str(''.join(list(filter(str.isalpha, minor))))+str(''.join(list(filter(str.isdigit, minor)))))
            
            Ymajor=SAGB_rereduced[SAGB_rereduced.Species.str.contains(major_search,case=False)]['ProdFact']
            Yminor=SAGB_rereduced[SAGB_rereduced.Species.str.contains(minor_search,case=False)]['ProdFact']
            
            if (Yminor.size>0) and (Ymajor.size>0):
                
                Ymajor=Ymajor.item()
                Yminor=Yminor.item()
                # d=(((ratio_solar*(10**Yminor)/(10**Ymajor)*dilu[i])/ratio_solar)-dilu[i])*1000
                # d=(((ratio_solar*(10**Yminor)/(10**Ymajor))/ratio_solar)-1)*1000
                d=((ratio_solar*(10**Yminor)/(10**Ymajor))/ratio_solar)

                # if np.isnan(d)==False:
                axs[h].scatter(x[i],d,s=1*j**3,marker='o',c=colors[i].reshape(1,-1), alpha=0.5, edgecolors='k')
                
                encountered_ratio.append(iso_ratio.ratio_name.iloc[i]) if iso_ratio.ratio_name.iloc[i] not in encountered_ratio else encountered_ratio
                tickposition.append(x[i]) if x[i] not in tickposition else tickposition
    
    
    
    
    AGB_reduced=AGB.loc[AGB.Metallicity==Z]
    
    for j in AGB_mass:
        AGB_rereduced=AGB_reduced.loc[AGB_reduced.Mass==j]
    
        for i in range(0,len(encountered_ratio)):
            extracted_rat=iso_ratio.loc[iso_ratio.ratio_name==encountered_ratio[i]]
            split_string = extracted_rat.ratio_name.item().split("/", 1)
            major=split_string[1]
            minor=split_string[0]
            ratio_solar=extracted_rat.ratio.item()
            
            major_search=str(str(''.join(list(filter(str.isalpha, major))))+str(''.join(list(filter(str.isdigit, major)))))
            minor_search=str(str(''.join(list(filter(str.isalpha, minor))))+str(''.join(list(filter(str.isdigit, minor)))))
    
                
            Ymajor=AGB_rereduced.filter(regex=str.lower(major_search))
            Yminor=AGB_rereduced.filter(regex=str.lower(minor_search))
            
            if (Yminor.size>0) and (Ymajor.size>0):
                Ymajor=Ymajor.values.item()
                Yminor=Yminor.values.item()
                # d=(((ratio_solar*(10**Yminor)/(10**Ymajor)*dilu[i])/ratio_solar)-dilu[i])*1000
                # d=(((Yminor/Ymajor)/ratio_solar)-1)*1000
                d=((Yminor/Ymajor)/ratio_solar)
             
                # if np.isnan(d)==False:
                axs[h].scatter(tickposition[i],d,s=1*j**3,marker='s',c=colors[i].reshape(1,-1), alpha=0.5, edgecolors='k')
                
    axs[h].set_title('Z= '+str(Z))
    axs[h].set_yscale('log')
    if h>1:
        # Set number of ticks for x-axis
        axs[h].set_xticks(tickposition)
        # ax5.axis([0,23,-500,25000])
        # Set ticks labels for x-axis
        axs[h].set_xticklabels(encountered_ratio, rotation='vertical', fontsize=18)



  