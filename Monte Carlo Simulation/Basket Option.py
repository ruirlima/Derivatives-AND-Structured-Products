import numpy as np
import pandas as pd
from pandas_datareader.data import DataReader
from datetime import date

#----------------------------------------------------ASSUMPTIONS--------------------------------------------------------
# Generic Assumptions
rate = 0.03                                 # discount rate
maturity = 1                                # Option maturity
days = 252                                  # days in a trading year
simulations = 1000                           # number of simulations
K = 1                                       # 1 for ATM option (=100%)

# Basket assumption
# Period to be used for covariance matrix
start_date = date(2020,1,1)
end_date = date(2020,7,31)
tickers = ['AMZN','FB','AAPL','GOOGL']
# This is multiplied by asset price, so DOES NOT represent asset percentage of overall basket.
# S1 * w1 + S2 * w2 + ..... + Sn * wn
# NOT
# sum(Si)*wi
weight = [0.5,0.2,0.2,0.1]

#-----------------------------------------------IMPORT PRICE DATA & RETURNS---------------------------------------------
# Get daily data from Yahoo
panel = DataReader(tickers,'yahoo', start_date, end_date)
S_close = panel['Adj Close'].unstack().unstack(level=0)
returns = pd.DataFrame(np.nan,columns=tickers,index=panel.index)
# Calculate stock returns
for i, idx in enumerate(returns.index):
    if i == 0:
        continue
    else:
        returns.iloc[i,:] = np.log(S_close.iloc[i,:] / S_close.iloc[i-1,:])

#-----------------------------------------------------Functions---------------------------------------------------------
def GBM(S,drift,vol,dt,rand):
    # rand - cdf standard normal
    # vol is already in "rand" due to cholesky decomposition
    return S * np.exp(drift*dt-0.5 * vol ** 2 * dt + np.sqrt(dt) * rand)

#---------------------------------------------CORRELATIONS / CHOLESKY---------------------------------------------------
# Variance Covariance matrix
VarCovar = np.multiply(np.cov(np.array(returns.iloc[1:,:]).astype(float),rowvar=False,ddof=0),days)        #remove first row because NaN
# Cholesky decomposition for correlation in MCS
cholesky = np.linalg.cholesky(VarCovar)
volatility = np.sqrt(np.diag(VarCovar))
#----------------------------------------------------Simulation---------------------------------------------------------
# Initialise list for asset prices
df_S = [[S_close.iloc[-1,i] for i in range(len(tickers))]*simulations for _ in range(days*maturity)]
# Convert list to dataframe
df_S = pd.DataFrame(np.array(df_S),columns=pd.MultiIndex.from_product([list(range(1,simulations+1)),tickers]),index=list(range(1,days*maturity+1)))

# Set seed for MCS
np.random.seed(1)
# Use geometric brownian motion to simulate prices across time.
# Random variable considers correlation among assets via Cholesky decomposition
# Cholesky includes asset volatility, so GBM function does not have vol in diffusion segment
for i,row in df_S.iterrows():
    if i == 1:
        continue
    else:
        # Loop trough simulations to simplify random variable analysis
        for sim in range(1,simulations+1):
            # random factor for basket option needs to consider correlation among assets, so multiply by cholesky
            rand = np.random.normal(0,1,len(tickers))
            rand_cholesky = np.matmul(rand,cholesky.transpose())
            df_S.loc[i,(sim,)] = GBM(np.array(df_S.loc[i-1,(sim,)]),rate,volatility,1/days,rand_cholesky)

print(df_S)

#----------------------------------------------Call and Put premium-----------------------------------------------------
# Compute index value at maturity
print(df_S.head())
print(df_S.tail())
final_index = [np.sum(np.multiply(np.array(df_S.loc[days,(sim,)])/np.array(df_S.loc[1,(sim,)]),weight)) for sim in range(1,simulations+1)]
# Option prices
call = np.exp(-rate*maturity)*sum(max(final_index[i]-K,0) for i in range(simulations))/simulations
put = np.exp(-rate*maturity)*sum(max(K-final_index[i],0) for i in range(simulations))/simulations

print('\t Call\t\tPut')
print('{0:>9s}{1:>10s}'.format(str(round(call,3)),
                                str(round(put,3))))