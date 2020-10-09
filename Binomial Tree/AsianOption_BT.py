import math as m
from itertools import product
def AsianOption_binomial(S,K,rate,div,t,T,sig,step,iopt):
    dt = (T - t) / step                                  # period per step
    erdt = m.exp(rate*(T-t))                                # compounding factor
    ermqdt = m.exp((rate - div) * dt)                    # dividend effect
    up = m.exp(sig * dt ** 0.5)                          # up multiplier
    down = 1 / up                                        # down multiplier
    prob_up = (ermqdt - down) / (up - down)              # up movement probability
    prob_down = 1 - prob_up                              # down probability
    # Calculate Stock prices for each step
    perm = list(product([1,0],repeat=step))              # possible paths for asset price
    share_avg = [0]*len(perm)                            # Average asset price for each path
    coefficient = [0]*len(perm)                          # binomial coefficient for each path (p^n)*(1-p)^(k-n)
    # Loop through each possible path
    j=0 # iterator for share_avg and coefficient
    for comb in perm:
        temp = [S]*(step+1)   #temporary list for stock price each path
        k=1 #iterator for temp
        # Lopp through each outcome inside path to determine whether stock up or down. Change temp accordingly.
        for i in comb:
            if i == 1:
                temp[k] = temp[k-1]*up
                k+=1
            else:
                temp[k] = temp[k-1]*down
                k+=1
        # Fill share_avg list with average stock price for path
        share_avg[j] = sum(temp)/len(temp)
        # Fill coefficient list withbinomial coefficient for path
        coefficient[j] = prob_up**sum(comb)*prob_down**(step-sum(comb))
        j+=1
    # Calculate Option Payoff at Maturity
    # Expected payoff is outcome times probability
    option_payoff = [max(iopt*(share_avg[i]-K),0)*coefficient[i] for i in range(len(share_avg))]
    # Present value of expected payoff
    pv_option = sum(option_payoff)/erdt
    return round(pv_option,3)

def main():
    # Validate Inputs
    while True:
        try:
            S = float(input('Stock price:= '))
            if S < 0.00:
                continue
            else:
                break
        except ValueError as e:
            print('Enter positive float!')
            continue
    while True:
        try:
            K = float(input('Strike Price:= '))
            if K < 0:
                continue
            else:
                break
        except ValueError as e:
            print('Enter positive float!')
            continue
    while True:
        try:
            rate = float(input('rate (decimal):= '))
            if rate < 0.00 or rate > 1.00:
                continue
            else:
                break
        except ValueError as e:
            print('Enter float between 0 and 1!')
            continue
    while True:
        try:
            div = float(input('dividend yield (decimal):= '))
            if div < 0.00 or div >1.00:
                continue
            else:
                break
        except ValueError as e:
            print('Enter float between 0 and 1!')
            continue
    while True:
        try:
            t = float(input('Time now (years):= '))
            if t < 0:
                continue
            else:
                break
        except ValueError as e:
            print('Enter positive float!')
            continue
    while True:
        try:
            T = float(input('Time maturity (years):= '))
            if T <= t:
                print('Time to maturity > time now')
                continue
            else:
                break
        except ValueError as e:
            print('Enter positive float larger than t')
            continue
    while True:
        try:
            sig = float(input('Standard deviation (decimal):= '))
            if sig < 0:
                continue
            else:
                break
        except ValueError as e:
            print('Enter positive float!')
            continue
    while True:
        try:
            step = int(input('Number of steps in binomial tree:= '))
            if step < 0:
                continue
            else:
                break
        except ValueError as e:
            print('Enter positive integer!')
            continue
    while True:
        try:
            iopt = int(input('iopt -> 1:= Call Option | -1:= Put Option: '))
            if iopt != 1 and iopt != -1:
                continue
            else:
                break
        except ValueError as e:
            print('1:= Call Option | -1:= Put Option:')
            continue

    # Determine which option type user chose, to put it in the output
    if iopt == 1:
        iopt_output = 'call'
    else:
        iopt_output = 'put'
    print('\nThe '+iopt_output+' option premium is: {0:>10s}'.format('£'+str(AsianOption_binomial(S,K,rate,div,t,T,sig,step,iopt))))

#--------------------------------------------------
main()