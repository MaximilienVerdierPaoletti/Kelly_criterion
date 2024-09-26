# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 15:18:49 2019

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


#Set working directory
os.chdir('D:/Work/Postdoc/Simulations/FRUITY/Dredge-ups/Extended 13Cpocket/')

DUlistfiles=os.listdir() # To retrieve the names of the files

for i in range(0,len(DUlistfiles)-1): #Loop over the number of files
    
    Mass=float(str(DUlistfiles[i][9]+'.'+DUlistfiles[i][11])) #retrieve mass out of the file name
    if DUlistfiles[i][13]=='s':#retrieve mass out of the file name
        Z=0.02 #Solar metallicity
    else:
        Z=int(DUlistfiles[i][13])*10**(-int(DUlistfiles[i][15]))
    
    
    X=pd.read_table(str(DUlistfiles[i]), header=None, delim_whitespace=True) #Do not recognize first line as header
    X=X.drop(X.index[0]) #Remove first line (bad header)
    X=X.dropna(how='all', axis=1) # Remove empty columns
    X.insert(0,"Mass",Mass) #Add the column mass
    X.insert(1,"Metallicity",Z) #Add the column metallicity
    if i==0: # If first iteration then create "data" variable
        data=X
    else: #if not, append new data to already existing variable
        data=data.append(X) 
        
        

head=['Mass','Metallicity','Isotope','A','Z','Initial','SDU'] # Beginning of header
TDU_num=len(data.columns)-len(head) #Retrieving the maximum number of TDU events in the dataset

for i in range(1,TDU_num+1): # Create the following of the header by adding the TDU-events
    head.append('TDU_'+str(i))
    
data.columns=head # Changing columns names

data.to_pickle('13CpockextDU_table.pkl') #save data as .pkl
data.to_excel('13CpockextDU_table.xlsx') #save data as .xlsx
