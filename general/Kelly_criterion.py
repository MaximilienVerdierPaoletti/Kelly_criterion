# -*- coding: utf-8 -*-

"""
Created on Thur 26 Sept 2024
@author: mverdier

"""

import numpy as np
import pandas as pd  # enables the use of dataframe
import matplotlib
import matplotlib.pyplot as plt
import math
from scipy.stats import binom


def close_factors(number):
    '''
    find the closest pair of factors for a given number
    '''
    factor1 = 0
    factor2 = number
    while factor1 + 1 <= factor2:
        factor1 += 1
        if number % factor1 == 0:
            factor2 = number // factor1

    return factor1, factor2

Wo=50 # Initial money
Ph=0.7 #Heads prob
Pt=0.3 # Tails prob
Payout=4/6 # Payout odds. Bet 6 get head win 4 else loose 6. Bet 3 get 2 else lose 3
N=100 # Number of flips
R=list(range(N+1)) # list of flips
N_bet=[2,5,6,10,50,100] # How much do you bet per flip
Nmax_games=40 # Number of games played


#------- Binomial distribution of the biased coin
dist = [binom.pmf(r, N, Ph) for r in R]
plt.figure()
plt.bar(R, dist)


# Wins and loss during each games
f1,f2=close_factors(len(N_bet))
fwin, axs = plt.subplots(f1,f2)
axs=axs.ravel()
fhist, axhist = plt.subplots(1,1)
for h in range(len(N_bet)):
    Wn_l = np.empty(0)
    Ngames = 1  # Initial game
    while Ngames <Nmax_games:
        HandT=np.random.binomial(n=1, p=Ph, size=100) # Get the Head and Tail distribution
        Wager=np.where(HandT==1, (1+Payout), -1) # Get the equivalent as the wager distribution
        Wn=np.array(Wo) # Initialized the array of wins and loss
        Wi=Wo
        for i in range(N):
            Wn=np.append(Wn,Wi+N_bet[h]*Wager[i])
            Wi=Wn[i+1]

        Wn_l=np.append(Wn_l,Wn[-1])

        axs[h].plot(R, Wn)
        axs[h].set_ylim(0,3000)
        axs[h].set_title('Bet of '+str(N_bet[h])+' €')
        Ngames += 1


    # ----- In kelly criterion : Conversion of the final wealth into a growth rate
    # Growth rate (g) : If Wo grew by g for each of N flips to yield Wn, what would g be ?
    # Wn = Wo(1+g)^N --> g = (Wn/Wo)^(1/N)-1
    xbins = np.linspace(-20,10,100)
    g=((Wn_l/Wo)**(1/N)-1)*100
    axhist.hist(g,bins = xbins, density = False,histtype = 'step',label=N_bet[h])

axhist.legend()
axhist.set_title('Growth rate (%)')
plt.show(block=True)





