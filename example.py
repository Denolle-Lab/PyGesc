from gesc import *
import os
os.system('clear')
## TIME /FREQUENCY:
typ = 't'          # 't': enters time series length (win= in s) and rate (rate=1/s)
                    # 'f': enters frequency series length (win=in Hz) and rate
                    # (1/Hz)
win = 500          # (in s) window length
rate = 0.05       # (in 1/s) sampling rate (10Hz)
FT = make_ft(win,rate,typ)
FT['Fmax']=2         # calculate until Fmax if < Nyquist
max_mode=1         # maximum number of surface-wave modes to ouput
freq=1/(1+np.arange(10)) # plots eigenfunctions for periods 1 -> 10s

## MEDIUM
H = np.array([0, 5, 20])         # (km) vector of interface depths 
beta = 3*np.array([0.5, 1])      # (km/s) S-wavespeed 
alpha = 5.4*np.array([0.5, 1])   # sqrt(3)*beta # (km/s) P-wavespeed 
rho = 2.27*np.array([0.5, 1])    # (kg/dm^3) density
Dmax=30*max(beta)/FT['df']       # (km) "half space" effective thickness, matters for long period
typ = ('cst','lin')              # constant within layers: 'cst', all gradients 'grad'
MED = make_layers(H,Dmax,alpha,beta,rho,typ)

NR = {'typ':'rec',       # choose recursive
      'N1':  30,            # lower bound for number of points / layer
      'N2': 50,
      'NR':  100}         # maximum number of points / layers

## GESC
bigC = gesc(MED,FT,NR,max_mode,freq) # plot eigenfunctions at specific  frequencies
