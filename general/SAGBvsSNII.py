# -*- coding: utf-8 -*-

"""
Created on Mon Sep 09 2024
@author: mverdier

"""

# try:
#     from IPython import get_ipython
#     get_ipython().run_line_magic('reset', '-f')
# except:
#     pass

#%% Modules


import os  # Import module that enable management of paths
import numpy as np
import pandas as pd  # enables the use of dataframe
# from scipy import sparse as sp
import matplotlib
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt  # Enables plotting of data
from mpl_toolkits import mplot3d
import pylab  # Enables the use of the 'savefig' function
# from random import randint
import glob  # enables to search for specific files in directory and save their names in a list
import re  # enables to remove part of strings and look for parts in it
# import sys
import math
import seaborn as sns
from nugridpy import nugridse as nse

sns.set_theme(style="whitegrid")

from scipy.stats import alpha


#%%

def mass_selSN(models_path):
    # Ask the user for the desired dataset
    print("\n" * 100)
    i = input('What dataset do you want to use: Heger&Woosley 2009 (0) or Rauscher et al. 2002 (1) or Woosley&Heger 2007 (2) or Pignatari 2015 (3)? \n')
    i = int(i)

    if i == 0:
        models = glob.glob(models_path[i] + 's*')
        models = [m.replace("\\", '/') for m in models]
        masses = models
        masses = [re.sub(models_path[i], '', m) for m in masses]
        masses = [re.sub('@presn', '', m) for m in masses]
        masses = [re.sub('s', '', m) for m in masses]

    elif i == 1:
        models = glob.glob(models_path[i] + 's*')
        models = [m.replace("\\", '/') for m in models]
        masses = models
        masses = [re.sub(models_path[i], '', m) for m in masses]
        masses = [re.sub('.presn_comp', '', m) for m in masses]
        masses = [re.sub('a28', '', m) for m in masses]
        masses = [re.sub('s', '', m) for m in masses]

    elif i == 2:
        models = glob.glob(models_path[i] + 's*')
        models = [m.replace("\\", '/') for m in models]
        masses = models
        masses = [re.sub(models_path[i], '', m) for m in masses]
        masses = [re.sub('.presn_comp', '', m) for m in masses]
        masses = [re.sub('a28', '', m) for m in masses]

    elif i == 3:
        models = glob.glob(models_path[i] + 'm*')
        models = [m.replace("\\", '/') for m in models]
        masses = models
        masses = [re.sub(models_path[i], '', m) for m in masses]
        masses = [re.sub('.presn_comp', '', m) for m in masses]
        masses = [re.sub('a28', '', m) for m in masses]

    # Ask the user for the desired simulation
    # print("\n"*100)
    if i !=3 :
        print("\n")
        print('Available model masses:', masses)
        mass = input('What mass for the star? \n')
        print("\n" * 5)
        # Retrieve the corresponding file name
        matching = [s for s in models if mass in s]
        if len(matching) > 1:
            print("More than one model for this mass : ")
            for (i, item) in enumerate(matching, start=0):
                print(item, ' select : ', i)
            idx = input('Which model to select ? \n')
            matching = matching[int(idx)]
    else: # For Pignatari models
        Pigna_idx = input('Which mass-loss model (type numer starting from 0): ' + str(masses) + '\n')
        Pigna_idx = int(Pigna_idx)
        matching = models[Pigna_idx]
        mass = 25 # Only 25 solar masses models

    return matching, i, mass


def extraction_cleaningSN(models_path):
    # Columns names for different simulations
    header_Heger = ['grid', 'celloutertotalmass', 'cellouterradius', 'celloutervelocity', 'celldensity',
                    'celltemperature', 'cellpressure', 'cellspecificenergy', 'cellspecificentropy', 'cellangularvelocity', 'cellA_bar',
                    'cellY_e', 'stability']
    header_woos = ['zone', 'mass below', 'zonemass']
    colnames_path = ['../nucleo_models/Heger08/hegercols.txt', '../nucleo_models/Rauscher2002/PreSN_composition/', '../nucleo_models/Woosley07/']

    matching, i, mass = mass_selSN(models_path)

    if i == 0:  # For the Heger 2009 data files
        # Heger files have column with spaces in their names. read_table always consider this as a delimiter and we have have therefore more names than columns. Leading to a misattribution of data
        data = pd.read_table(matching[0], sep=r"\s+", header=1)  # Open the file once reading automatically the names of column
        colnames = data.columns.tolist()  # Save the names of columns in a list
        ind = data.columns.get_loc("nt1")  # Look for the index of the column named "nt1", columns prior that one contains space that's problematic for determining the real number of columns when read
        del colnames[0:ind]  # Delete all columns which names contain spaces
        [colnames.insert(j, header_Heger[j]) for j in range(0, len(header_Heger))]  # Inserts the rights names of the first columns
        del data, ind, header_Heger  # delete all obsolete variables
        data = pd.read_table(matching[0], names=colnames, sep=r'\s+', index_col=False, skiprows=[0, 1])  # remove the two first lines that are wrong headers

    elif i == 1:  # Rauscher 2002
        colnames = open(colnames_path[i] + mass + '_colnames.txt').read().split()
        data = pd.read_table(matching[0], sep=r"\s+", names=colnames, index_col=False)  # index_col=False enable to automatically index rows rather than giving them a column gathering their names
        data = data.drop([0], axis=0)  # Remove only the first line
        data = data.drop(data.tail(1).index, axis=0)  # remove last line because mass coordinate is close to core
        data.columns = data.columns.str.replace(' ', '')  # remove space from column names

    elif i == 2:  # Woosley&Heger 2007
        data = pd.read_table(matching[0], sep=r"\s+", header=0)  # Open the file once reading automatically the names of column
        colnames = data.columns.tolist()  # Save the names of columns in a list
        ind = data.columns.get_loc("nt1")  # Look for the index of the column named "nt1", that is the first column after all the ones which names contain spaces
        del colnames[0:ind]  # Delete all columns which names contain spaces
        [colnames.insert(j, header_woos[j]) for j in range(0, len(header_woos))]  # Inserts the rights names of the first columns
        del data, ind, header_woos  # delete all obsolete variables
        data = pd.read_table(matching[0], sep=r"\s+", header=1, names=colnames)
        data.columns = data.columns.str.replace(' ', '')  # remove space from column names

    elif i == 3:
        data = pd.read_table(str(matching), sep=" ", header=0)
        data.columns = data.columns.str.lower().str.replace('_', '')
        data.columns = data.columns.str.replace('mass', 'Mass')

    # Convert all the dataframe from object to floats
    cols = data.columns[data.dtypes.eq(object)]
    data[cols] = data[cols].apply(pd.to_numeric, errors='coerce', axis=0)  # /!\ this is what takes so long

    # Conversion of mass from gram to solar masses
    if i == 0:
        M = data.celloutertotalmass / 2E33  # Divided by the mass of the sun in grams as units are grams in the dataset
        data = pd.concat([pd.DataFrame({'Mass': M}),data],axis=1)
    elif i in [1,2]:
        M = data.massbelow / 2E33
        data = pd.concat([pd.DataFrame({'Mass': M}), data], axis=1)
        
    return data, i


def zone_definerSN(data, j):
    zones=[]
    if j !=3: # If not Pignatari models
        interest=['h1','he4','c12','n14','o16','ne20','si28','s32']
    else:
        interest = ['c12', 'n14', 'o16', 'si28', 's32']
    data_sel=data[interest]
    for i in data_sel.index:
        largest_elem=data_sel.loc[i].sort_values().index[-2:]
        if all(x in largest_elem for x in ['s32', 'si28']): zones = zones + ['Si/S']
        elif all(x in largest_elem for x in ['si28', 'o16']): zones = zones +['O/Si']
        elif all(x in largest_elem for x in ['si28', 'c12']):zones = zones + ['C/Si'] # Only in Pignatari models
        elif all(x in largest_elem for x in ['o16', 'ne20']) : zones = zones + ['O/Ne']
        elif all(x in largest_elem for x in ['o16', 'c12']): zones = zones +['C/O']
        elif all(x in largest_elem for x in ['o16', 'he4']) : zones = zones + ['O/He']
        elif all(x in largest_elem for x in ['c12', 'he4']): zones = zones + ['He/C']
        elif all(x in largest_elem for x in ['n14', 'he4']): zones = zones + ['He/N']
        elif data_sel.loc[i].sort_values().index[-1] == "h1": zones = zones + ['H/He']
        else: zones+ ['undefined']

        Z=pd.DataFrame({'Regions': zones})

    return pd.concat([Z,data],axis=1)


def SAGB_model(model_path):
    models_names = ['VW93', 'VW-M', 'testcases']
    # Extract models data
    pd.set_option("display.precision", 8)  # Fix the number of significant figures to 8
    Mini = pd.read_excel(model_path + 'InitialCompositions.xlsx', header=0)  # Extract initial compositions of the stars model
    m = input('Which mass-loss model: ' + str(models_names) + '(0,1,2)\n')  # Selectiom by the user of the desired model
    m = int(m)  # Conversion of answer into integer
    mod = pd.read_excel(model_path + '/reformatted' + str(models_names[m]) + '.xlsx', header=0)  # "E-" columns correspond to yields corresponding to extrapolated Thermal Pulses (TP) after
    # convergence of the main model
    starmass = mod.Mass.unique()  # Retrieve the different star mass available
    starmet = mod.Metallicity.unique()  # Retrieve the different metallicities available

    return Mini, mod, starmass, starmet

#%%

if __name__ == "__main__":

    # ----------------------------- Supernovae Models  -------------------------------------------#
    models_pathSN = ['../nucleo_models/Heger08/preSN/', '../nucleo_models/Rauscher2002/PreSN_composition/',
                     '../nucleo_models/Woosley07/preSN_composition/s*','../nucleo_models/Pignatari-models/',]
    model_pathSAGB = '../nucleo_models/Doherty15(Super-AGB)/'

    #Solar ratios (Lodders 2003)
    RS1312C=0.011237
    RS1716O=0.00037305
    RS1816O=0.002004965
    RS2524Mg=0.126597989
    RS2624Mg=0.139381904
    RS2928Si=0.050775236
    RS3028Si=0.033470671
    RS5456Fe=0.063702945
    RS5756Fe=0.023094361
    RS6058Ni=0.38519821
    RS6158Ni=0.016744299
    RS6258Ni=0.053388154

    # ------------- Supernovae

    dataSN, ind = extraction_cleaningSN(models_pathSN)
    dataSN = zone_definerSN(dataSN,ind)

    fig, [axSN,axReg]= plt.subplots(1,2,figsize=(12, 12))
    # fig, axSN = plt.subplots(1,1)
    lw = 2
    if ind !=3:
        axSN.semilogy(dataSN.Mass, dataSN.h1, label='1H', linewidth=lw)
        axSN.semilogy(dataSN.Mass, dataSN.he4, label='4He', linewidth=lw)
        axSN.semilogy(dataSN.Mass, dataSN.ne20, label='20Ne', linewidth=lw)
    axSN.semilogy(dataSN.Mass, dataSN.c12, label='12C', linewidth=lw)
    axSN.semilogy(dataSN.Mass, dataSN.n14, label='14N', linewidth=lw)
    axSN.semilogy(dataSN.Mass, dataSN.o16, label='16O', linewidth=lw)
    axSN.semilogy(dataSN.Mass, dataSN.si28, label='28Si', linewidth=lw)
    axSN.semilogy(dataSN.Mass, dataSN.s32, label='32S', linewidth=lw)

    axSN.set_xlim(2, 12)
    axSN.set_ylim(10 ** -8, 1)
    axSN.legend()

    axReg.semilogy(dataSN.Mass, dataSN.fe54, label='54Fe', linewidth=lw)
    axReg.semilogy(dataSN.Mass, dataSN.fe56, label='56Fe', linewidth=lw)
    axReg.semilogy(dataSN.Mass, dataSN.ni58, label='58Ni', linewidth=lw)
    axReg.semilogy(dataSN.Mass, dataSN.ni60, label='60Ni', linewidth=lw)
    axReg.semilogy(dataSN.Mass, dataSN.ni61, label='61Ni', linewidth=lw)
    axReg.semilogy(dataSN.Mass, dataSN.ni62, label='62Ni', linewidth=lw)
    axReg.set_xlim(2, 12)
    axReg.set_ylim(10 ** -8, 1)
    axReg.legend()

    # if ind !=3:
    #     g = sns.lineplot(data=dataSN, x="Mass", y="h1", hue="Regions",ax=axReg)
    #     g = sns.lineplot(data=dataSN, x="Mass", y="he4", hue="Regions",ax=axReg)
    # g = sns.lineplot(data=dataSN, x="Mass", y="si28", hue="Regions",ax=axReg)
    # g = sns.lineplot(data=dataSN, x="Mass", y="o16", hue="Regions",ax=axReg)
    # g.set(yscale='log')
    # axReg.set_xlim(2, 12)
    # axReg.set_ylim(10 ** -8, 1)

    dataSN['R1716O_SN'] = ((dataSN.o17/dataSN.o16)/RS1716O-1)*1000
    dataSN['R1816O_SN'] = ((dataSN.o18 / dataSN.o16) / RS1816O - 1) * 1000
    dataSN['R2524Mg_SN'] = ((dataSN.mg25 / dataSN.mg24) / RS2524Mg - 1) * 1000
    dataSN['R2624Mg_SN'] = ((dataSN.mg26 / dataSN.mg24) / RS2624Mg - 1) * 1000
    dataSN['R5456Fe_SN'] = ((dataSN.fe54 / dataSN.fe56) / RS5456Fe - 1) * 1000
    dataSN['R6058Ni_SN'] = ((dataSN.ni60 / dataSN.ni58) / RS6058Ni - 1) * 1000
    dataSN['R6158Ni_SN'] = ((dataSN.ni61 / dataSN.ni58) / RS6158Ni - 1) * 1000
    dataSN['R6258Ni_SN'] = ((dataSN.ni62 / dataSN.ni58) / RS6258Ni - 1) * 1000

    f_deltaSN = plt.figure()
    plt.semilogy(dataSN.Mass, dataSN.R1716O_SN, label='d-17O/16O', linewidth=lw)
    plt.semilogy(dataSN.Mass, dataSN.R1816O_SN, label='d-18O/16O', linewidth=lw)
    plt.semilogy(dataSN.Mass, dataSN.R2524Mg_SN, label='d-25Mg/24Mg', linewidth=lw)
    plt.semilogy(dataSN.Mass, dataSN.R2624Mg_SN, label='d-26Mg/24Mg', linewidth=lw)
    plt.semilogy(dataSN.Mass, dataSN.R5456Fe_SN, label='d-54Fe/56Fe', linewidth=lw)
    plt.semilogy(dataSN.Mass, dataSN.R6058Ni_SN, label='d-60Ni/58Ni', linewidth=lw)
    plt.semilogy(dataSN.Mass, dataSN.R6158Ni_SN, label='d-61Ni/58Ni', linewidth=lw)
    plt.semilogy(dataSN.Mass, dataSN.R6258Ni_SN, label='d-62Ni/58Ni', linewidth=lw)
    # plt.axis([2, 12,10 ** -8, 1])
    plt.legend()

    CSi = dataSN.loc[dataSN.Mass.between(6.8145, 6.8247, inclusive='both')]
    Onova = dataSN.loc[dataSN.Mass.between(7.05, 7.2, inclusive='both')]
    
    # ------------- S-AGB

    Mini, mod, starmass, starmet = SAGB_model(model_pathSAGB)

    # R1716O = {'Yield': mod.loc[mod.Species=='o17']['MassExp'].values/mod.loc[mod.Species=='o16']['MassExp'].values, 'Species': '17O'}
    # R1816O = {'Yield': mod.loc[mod.Species=='o18']['MassExp'].values/mod.loc[mod.Species=='o16']['MassExp'].values, 'Species': '18O'}
    # R2524Mg = {'Yield': mod.loc[mod.Species=='mg25']['MassExp'].values/mod.loc[mod.Species=='mg24']['MassExp'].values, 'Species': '25Mg'}
    # R2624Mg = {'Yield': mod.loc[mod.Species=='mg26']['MassExp'].values/mod.loc[mod.Species=='mg24']['MassExp'].values, 'Species': '26Mg'}
    # R5456Fe = {'Yield': mod.loc[mod.Species=='fe54']['MassExp'].values/mod.loc[mod.Species=='fe58']['MassExp'].values, 'Species': '54Fe'}
    # R5756Fe = {'Yield': mod.loc[mod.Species=='fe57']['MassExp'].values/mod.loc[mod.Species=='fe58']['MassExp'].values, 'Species': '57Fe'}
    # R6058Ni = {'Yield': mod.loc[mod.Species=='ni60']['MassExp'].values/mod.loc[mod.Species=='ni58']['MassExp'].values, 'Species': '60Ni'}
    # R6158Ni = {'Yield': mod.loc[mod.Species=='ni61']['MassExp'].values/mod.loc[mod.Species=='ni58']['MassExp'].values, 'Species': '61Ni'}
    # R6258Ni = {'Yield': mod.loc[mod.Species=='ni62']['MassExp'].values/mod.loc[mod.Species=='ni58']['MassExp'].values, 'Species': '62Ni'}

    # data_SAGB=pd.concat([mod.loc[mod.Species=='ni62',['Mass', 'Metallicity']].reset_index(),
    #                      R1716O, R1816O, R2524Mg,R2624Mg, R5456Fe, R5756Fe, R6058Ni, R6158Ni, R6258Ni ],axis=1)

    R1716O = {'Yield': ((mod.loc[mod.Species=='o17']['MassExp'].values/mod.loc[mod.Species=='o16']['MassExp'].values)/RS1716O-1)*1000  , 'Species': '17O'}
    R1816O = {'Yield': ((mod.loc[mod.Species=='o18']['MassExp'].values/mod.loc[mod.Species=='o16']['MassExp'].values)/RS1816O-1)*1000, 'Species': '18O'}
    R2524Mg = {'Yield': ((mod.loc[mod.Species=='mg25']['MassExp'].values/mod.loc[mod.Species=='mg24']['MassExp'].values)/RS2524Mg-1)*1000, 'Species': '25Mg'}
    R2624Mg = {'Yield': ((mod.loc[mod.Species=='mg26']['MassExp'].values/mod.loc[mod.Species=='mg24']['MassExp'].values)/RS2624Mg-1)*1000, 'Species': '26Mg'}
    R5456Fe = {'Yield': ((mod.loc[mod.Species=='fe54']['MassExp'].values/mod.loc[mod.Species=='fe58']['MassExp'].values)/RS5456Fe-1)*1000, 'Species': '54Fe'}
    R5756Fe = {'Yield': ((mod.loc[mod.Species=='fe57']['MassExp'].values/mod.loc[mod.Species=='fe58']['MassExp'].values)/RS5756Fe-1)*1000, 'Species': '57Fe'}
    R6058Ni = {'Yield': ((mod.loc[mod.Species=='ni60']['MassExp'].values/mod.loc[mod.Species=='ni58']['MassExp'].values)/RS6058Ni-1)*1000, 'Species': '60Ni'}
    R6158Ni = {'Yield': ((mod.loc[mod.Species=='ni61']['MassExp'].values/mod.loc[mod.Species=='ni58']['MassExp'].values)/RS6158Ni-1)*1000, 'Species': '61Ni'}
    R6258Ni = {'Yield': ((mod.loc[mod.Species=='ni62']['MassExp'].values/mod.loc[mod.Species=='ni58']['MassExp'].values)/RS6258Ni-1)*1000, 'Species': '62Ni'}

    data_SAGB = pd.concat([pd.concat([mod.loc[mod.Species == 'o17', ['Mass', 'Metallicity']].reset_index(), pd.DataFrame(R1716O)], axis=1),
               pd.concat([mod.loc[mod.Species == 'o18', ['Mass', 'Metallicity']].reset_index(), pd.DataFrame(R1816O)], axis=1),
               pd.concat([mod.loc[mod.Species == 'mg25', ['Mass', 'Metallicity']].reset_index(), pd.DataFrame(R2524Mg)], axis=1),
               pd.concat([mod.loc[mod.Species == 'mg26', ['Mass', 'Metallicity']].reset_index(), pd.DataFrame(R2624Mg)], axis=1),
               pd.concat([mod.loc[mod.Species == 'fe54', ['Mass', 'Metallicity']].reset_index(), pd.DataFrame(R5456Fe)], axis=1),
               pd.concat([mod.loc[mod.Species == 'fe57', ['Mass', 'Metallicity']].reset_index(), pd.DataFrame(R5756Fe)], axis=1),
               pd.concat([mod.loc[mod.Species == 'ni60', ['Mass', 'Metallicity']].reset_index(), pd.DataFrame(R6058Ni)], axis=1),
               pd.concat([mod.loc[mod.Species == 'ni61', ['Mass', 'Metallicity']].reset_index(), pd.DataFrame(R6158Ni)], axis=1),
               pd.concat([mod.loc[mod.Species == 'ni62', ['Mass', 'Metallicity']].reset_index(), pd.DataFrame(R6258Ni)], axis=1)],axis=0)


    plt.figure()
    g=sns.scatterplot(
        data=data_SAGB, x="Species", y="Yield", hue="Metallicity", size="Mass",
        sizes=(30, 300), legend="full", palette="flare")
    g.set(yscale='symlog')
    plt.ylabel('Delta values', fontweight='bold')
    plt.title('Super AGB Yields',fontsize=20)


    ExR1716O = {'Yield': ((mod.loc[mod.Species=='o17']['E-MassExp'].values/mod.loc[mod.Species=='o16']['E-MassExp'].values)/RS1716O-1)*1000  , 'Species': '17O'}
    ExR1816O = {'Yield': ((mod.loc[mod.Species=='o18']['E-MassExp'].values/mod.loc[mod.Species=='o16']['E-MassExp'].values)/RS1816O-1)*1000, 'Species': '18O'}
    ExR2524Mg = {'Yield': ((mod.loc[mod.Species=='mg25']['E-MassExp'].values/mod.loc[mod.Species=='mg24']['E-MassExp'].values)/RS2524Mg-1)*1000, 'Species': '25Mg'}
    ExR2624Mg = {'Yield': ((mod.loc[mod.Species=='mg26']['E-MassExp'].values/mod.loc[mod.Species=='mg24']['E-MassExp'].values)/RS2624Mg-1)*1000, 'Species': '26Mg'}
    ExR5456Fe = {'Yield': ((mod.loc[mod.Species=='fe54']['E-MassExp'].values/mod.loc[mod.Species=='fe58']['E-MassExp'].values)/RS5456Fe-1)*1000, 'Species': '54Fe'}
    ExR5756Fe = {'Yield': ((mod.loc[mod.Species=='fe57']['E-MassExp'].values/mod.loc[mod.Species=='fe58']['E-MassExp'].values)/RS5756Fe-1)*1000, 'Species': '57Fe'}
    ExR6058Ni = {'Yield': ((mod.loc[mod.Species=='ni60']['E-MassExp'].values/mod.loc[mod.Species=='ni58']['E-MassExp'].values)/RS6058Ni-1)*1000, 'Species': '60Ni'}
    ExR6158Ni = {'Yield': ((mod.loc[mod.Species=='ni61']['E-MassExp'].values/mod.loc[mod.Species=='ni58']['E-MassExp'].values)/RS6158Ni-1)*1000, 'Species': '61Ni'}
    ExR6258Ni = {'Yield': ((mod.loc[mod.Species=='ni62']['E-MassExp'].values/mod.loc[mod.Species=='ni58']['E-MassExp'].values)/RS6258Ni-1)*1000, 'Species': '62Ni'}

    data_ExSAGB = pd.concat([pd.concat([mod.loc[mod.Species == 'o17', ['Mass', 'Metallicity']].reset_index(), pd.DataFrame(R1716O)], axis=1),
               pd.concat([mod.loc[mod.Species == 'o18', ['Mass', 'Metallicity']].reset_index(), pd.DataFrame(R1816O)], axis=1),
               pd.concat([mod.loc[mod.Species == 'mg25', ['Mass', 'Metallicity']].reset_index(), pd.DataFrame(R2524Mg)], axis=1),
               pd.concat([mod.loc[mod.Species == 'mg26', ['Mass', 'Metallicity']].reset_index(), pd.DataFrame(R2624Mg)], axis=1),
               pd.concat([mod.loc[mod.Species == 'fe54', ['Mass', 'Metallicity']].reset_index(), pd.DataFrame(R5456Fe)], axis=1),
               pd.concat([mod.loc[mod.Species == 'fe57', ['Mass', 'Metallicity']].reset_index(), pd.DataFrame(R5756Fe)], axis=1),
               pd.concat([mod.loc[mod.Species == 'ni60', ['Mass', 'Metallicity']].reset_index(), pd.DataFrame(R6058Ni)], axis=1),
               pd.concat([mod.loc[mod.Species == 'ni61', ['Mass', 'Metallicity']].reset_index(), pd.DataFrame(R6158Ni)], axis=1),
               pd.concat([mod.loc[mod.Species == 'ni62', ['Mass', 'Metallicity']].reset_index(), pd.DataFrame(R6258Ni)], axis=1)],axis=0)

    plt.figure()
    gex=sns.scatterplot(
        data=data_ExSAGB, x="Species", y="Yield", hue="Metallicity", size="Mass",
        sizes=(30, 300), legend="full", palette="flare")
    gex.set(yscale='symlog')
    plt.ylabel('Delta values', fontweight='bold')
    plt.title('Super AGB Yields Extrapolated TP',fontsize=20)


    # plt.figure()
    # ax = plt.axes(projection="3d")
    # ax.plot_trisurf(mod.loc[mod.Species == 'o17']['Mass'],mod.loc[mod.Species == 'o17']['Metallicity'],(R17/RS1716O-1)*1000,
    #                 linewidth=0, antialiased=True,alpha=0.7)
    # ax.plot_trisurf(mod.loc[mod.Species == 'o17']['Mass'], mod.loc[mod.Species == 'o17']['Metallicity'], (R18/RS1816O-1)*1000, linewidth=0, antialiased=False)
    #
    # # Adding labels
    # ax.set_xlabel('Mass', fontweight='bold')
    # ax.set_ylabel('Yield', fontweight='bold')
    # ax.set_zlabel('Ratio Yield', fontweight='bold')
    plt.show(block=False)
    plt.show()


# %%
