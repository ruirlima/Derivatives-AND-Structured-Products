'''
'Computation Barrier Option Premium via Monte Carlo Simulation
'   Up-and-In
'   Up-and-out
'   Down-and-in
'   Down-and-out
'''
# Import Libraries
import pandas as pd
import numpy as np

#--------------------------------------------Assumptions/Variables------------------------------------------------------
S = 100                              #Initial underlying price
K = 100                              #Strike price
L = 110                              #Barrier level
rate = 0.03                          #risk free rate (discount rate)
vol = 0.2                            #volatility underlying
time = 1                             #Time to expiry in years
simulations = 5000                   #Number simulations in MCS
days_per_year = 252                  #Number of days per year

#-----------------------------------------------------Functions---------------------------------------------------------
def GBM(S,drift,vol,dt,rand):
    # rand - pdf standard normal
    return S * np.exp(drift*dt-0.5 * vol ** 2 * dt + vol * np.sqrt(dt) * rand)
#--------------------------------------------------Price Simulation-----------------------------------------------------
# Initialise dataframe for Stock prices
df_S = [[S]*simulations for _ in range(days_per_year*time)]
df_S = pd.DataFrame(df_S,index=[i for i in range(1,days_per_year*time+1)],columns=[i for i in range(1,simulations+1)])
# Set seed for cdf standard normal
np.random.seed(1)
# Use geometric brownian motion to simulate prices across tim
for i,row in df_S.iterrows():
    if i == 1:
        continue
    rand = np.random.normal(0,1,simulations)
    df_S.loc[i,:] = GBM(df_S.loc[i-1,:],rate,vol,1/days_per_year,rand)
print(df_S)
#--------------------------------------------------Option Payoff--------------------------------------------------------
# Call premium for different call possibilities
call_up_and_in = np.exp(-rate*time)*sum(max(df_S.iloc[-1,i]-K,0) for i in range(simulations) if df_S.iloc[:,i].max() > L)/simulations
call_up_and_out = np.exp(-rate*time)*sum(max(df_S.iloc[-1,i]-K,0) for i in range(simulations) if df_S.iloc[:,i].max() < L)/simulations
call_down_and_in = np.exp(-rate*time)*sum(max(df_S.iloc[-1,i]-K,0) for i in range(simulations) if df_S.iloc[:,i].min() < L)/simulations
call_down_and_out = np.exp(-rate*time)*sum(max(df_S.iloc[-1,i]-K,0) for i in range(simulations) if df_S.iloc[:,i].min() > L)/simulations
# Put premium for different call possibilities
put_up_and_in = np.exp(-rate*time)*sum(max(K-df_S.iloc[-1,i],0) for i in range(simulations) if df_S.iloc[:,i].max() > L)/simulations
put_up_and_out = np.exp(-rate*time)*sum(max(K-df_S.iloc[-1,i],0) for i in range(simulations) if df_S.iloc[:,i].max() < L)/simulations
put_down_and_in = np.exp(-rate*time)*sum(max(K-df_S.iloc[-1,i],0) for i in range(simulations) if df_S.iloc[:,i].min() < L)/simulations
put_down_and_out = np.exp(-rate*time)*sum(max(K-df_S.iloc[-1,i],0) for i in range(simulations) if df_S.iloc[:,i].min() > L)/simulations

#--------------------------------------------------Print--------------------------------------------------------
print('\n\t\t Up-and-In\t\tUp-and-Out\t\tDown-and-In\t\tDown-and-Out')
print('Call{0:>14s}{1:>16s}{2:>17s}{3:>17s}'.format(str(round(call_up_and_in,3)),
                                                str(round(call_up_and_out,3)),
                                                str(round(call_down_and_in,3)),
                                                str(round(call_down_and_out,3))))
print('Put{0:>15s}{1:>16s}{2:>17s}{3:>17s}\n'.format(str(round(put_up_and_in,3)),
                                                str(round(put_up_and_out,3)),
                                                str(round(put_down_and_in,3)),
                                                str(round(put_down_and_out,3))))


