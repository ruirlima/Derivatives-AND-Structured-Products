import numpy as np
import pandas as pd
import datetime as dt
import math as m
from pandas_datareader.data import DataReader
from scipy.stats import norm
import seaborn as sns
import matplotlib.pyplot as plt

#-----------------------------------------------FUNCTIONS---------------------------------------------
def GBM(S,drift,vol,dt,rand):
    # rand - cdf standard normal
    # vol is already in "rand" due to cholesky decomposition
    return S * np.exp(drift*dt-0.5 * vol ** 2 * dt + np.sqrt(dt) * rand)

def BlackScholes(S,K,vol,rate,dt,iopt):
    d1 = (np.log(S/K)+(rate+0.5*np.power(vol,2))*dt)/(vol*np.sqrt(dt))
    d2 = d1 - vol*np.sqrt(dt)
    premium = iopt*S*norm.cdf(iopt*d1) - iopt*K*np.exp(-rate*dt)*norm.cdf(iopt*d2)
    return premium

def Greeks(S,K,vol,rate,dt,iopt):
    d1 = (np.log(S/K)+(rate+0.5*np.power(vol,2))*dt)/(vol*np.sqrt(dt))
    delta = iopt*norm.cdf(iopt*d1)
    gamma = norm.pdf(d1)/(vol*S*np.sqrt(dt))
    return [delta,gamma]

#----------------------------------------------------ASSUMPTIONS--------------------------------------------------------
# Generic Assumptions
rate = 0.03                                 # discount rate
days = 252                                  # days in a trading year
simulations = 1000                          # number of simulations

# Options assumption
# Period to be used for covariance matrix
start_date = dt.date(2014,1,2)              # period for correlation estimation (start)
end_date = dt.date(2019,11,2)               # period for correlation estimation (end)
tickers = ['HLT','MAR']                     # underlying asset names in portfolio
K = {'HLT':0,'MAR':0}                                   # 0 for ATM. Otherwise insert value
time_maturity = 0.5                                     # time to maturity in years
notional = {'HLT':5000,'MAR':5000}                      # notional in options per asset
type = {'HLT':1,'MAR':-1}                               # 1:call    -1:put

# VaR
VaR_period = 25                             # period in days for Value at Risk
percentile = 5                              #

#------------------------------------------IMPORT PRICE DATA & RETURNS--------------------------------------------------
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
#---------------------------------------------CORRELATIONS / CHOLESKY---------------------------------------------------
# Variance Covariance matrix
VarCovar = np.multiply(np.cov(np.array(returns.iloc[1:,:]).astype(float),rowvar=False,ddof=0),days)        #remove first row because NaN
# Cholesky decomposition for correlation in MCS
cholesky = np.linalg.cholesky(VarCovar)
volatility = np.sqrt(np.diag(VarCovar))
rate_return = np.power(returns.mean(skipna=True)+1,days)-1

#-------------------------------------------------Simulation------------------------------------------------------------
# Initialise list for asset prices
df_S = [[S_close.iloc[-1,i] for i in range(len(tickers))]*simulations for _ in range(2)]
# Convert list to dataframe
df_S = pd.DataFrame(np.array(df_S),columns=pd.MultiIndex.from_product([list(range(1,simulations+1)),tickers]),index=list(range(1,3)))

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
            df_S.loc[i,(sim,)] = GBM(np.array(df_S.loc[i-1,(sim,)]),rate,volatility,VaR_period/days,rand_cholesky)

#-------------------------------------------------Option Premium--------------------------------------------------------
# Set strike price for ATM case
for key,value in K.items():
    if value == 0:
        K[key] = df_S.loc[1,(1,key)]

# Initiate option premium dataframe
option_BS = df_S.copy()
# Maturities for Black Scholes. time 0, and time VaR
maturity = [time_maturity,time_maturity-VaR_period/days]
for i,row in option_BS.iterrows():
    # Loop trough simulations to simplify random variable analysis
    for sim in range(1,simulations+1):
        option_BS.loc[i,(sim,)] = BlackScholes(S=np.array(option_BS.loc[i,(sim,)]),
                                               K=np.array(list(K.values())),
                                               vol = volatility,
                                               rate = rate,
                                               dt = maturity[i-1],
                                               iopt = np.array(list(type.values())))
# Change premium via Black Scholes
dC_BS = option_BS.loc[2,]-option_BS.loc[1,]

# Calculate greeks for delta and delta gamma approximation
greek = Greeks(S=np.array(df_S.loc[1,(1,)]),
               K=np.array(list(K.values())),
               vol=volatility,
               rate=rate,
               dt=time_maturity,
               iopt=np.array(list(type.values())))

# Change asset price
deltaS = df_S.loc[2,]-df_S.loc[1,]

# Delta approximation
dC_delta = dC_BS.copy()                          # initialise delta approximation dataframe
for sim in range(1,simulations+1):
    dC_delta.loc[sim,] = np.multiply(np.array(deltaS.loc[sim,]),greek[0])

# Delta-Gamma approximation
dC_delta_gamma = dC_delta.copy()                 # initialise delta gamma approximation dataframe
for sim in range(1,simulations+1):
    dC_delta_gamma.loc[sim,] = np.multiply(np.array(deltaS.loc[sim,]),greek[0]) + \
                               0.5*np.multiply(np.power(np.array(deltaS.loc[sim,]),2),greek[1])

#-------------------------------------------------Value at Risk---------------------------------------------------------
# Determine how many units of each asset
units = notional.copy()
for key,value in notional.items():
    units[key]=value/option_BS.loc[1,(1,key)]

# Determine distribution portfolio value in each method
dC_distribution = pd.DataFrame(0,columns=list(range(1,simulations+1)),index=['Black Scholes','Delta','Delta Gamma'])
for i in range(1,simulations+1):
    dC_distribution.loc['Black Scholes',i] = (dC_BS.loc[i,]*np.array(list(units.values()))).sum()
    dC_distribution.loc['Delta', i] = (dC_delta.loc[i,] * np.array(list(units.values()))).sum()
    dC_distribution.loc['Delta Gamma', i] = (dC_delta_gamma.loc[i,] * np.array(list(units.values()))).sum()

# Calculate VaR at certain percentile
VaR = pd.DataFrame(0,columns=['VaR'],index=['Black Scholes','Delta','Delta Gamma'])
VaR.loc['Black Scholes','VaR'] = np.percentile(dC_distribution.loc['Black Scholes',],percentile)
VaR.loc['Delta','VaR'] = np.percentile(dC_distribution.loc['Delta',],percentile)
VaR.loc['Delta Gamma','VaR'] = np.percentile(dC_distribution.loc['Delta Gamma',],percentile)

print(VaR)

#-------------------------------------------------PLOT---------------------------------------------------------
plt.hist(dC_distribution.loc['Black Scholes',],bins=round(m.sqrt(simulations)),color='r',alpha=0.1)
plt.hist(dC_distribution.loc['Delta',],bins=round(m.sqrt(simulations)),color='b',alpha=0.1)
plt.hist(dC_distribution.loc['Delta Gamma',],bins=round(m.sqrt(simulations)),color='g',alpha=0.1)
plt.show()