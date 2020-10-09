import math as m
def option_binomial(S,X,rate,div,t,T,sig,step,iopt,type):
    dt = (T - t) / step                                  # period per step
    erdt = m.exp(rate*dt)                                # compounding factor
    ermqdt = m.exp((rate - div) * dt)                    # dividend effect
    up = m.exp(sig * dt ** 0.5)                          # up multiplier
    down = 1 / up                                        # down multiplier
    prob_up = (ermqdt - down) / (up - down)              # up movement probability
    prob_down = 1 - prob_up                              # down probability
    size = step + 1                                      # Table size
    # Calculate Stock prices for each step
    share_table = [[0]*size for _ in range(size)]        #list comprehension - nested list to initialize share price matrix
    for i in range(0,size):
        for j in range(0,i+1):
            factor = m.pow(up,j)*m.pow(down,i-j)
            share_table[i][j]= round(factor * S,7)

    # Calculate Option Payoff at Maturity
    option_payoff = [0]*size
    for i in range(0,size):
        option_payoff[i]=round(max(iopt*(share_table[-1][i]-X),0),4)

    # Calculate Option Payoff each Step
    option_table = [[0]*size for _ in range(size)]   #list comprehension - nested list to initialize share price matrix
    option_table[step] = option_payoff
    for i in range(step-1,-1,-1):
        for j in range(i,-1,-1):
            option_table[i][j] = round((prob_up*option_table[i+1][j+1]+prob_down*option_table[i+1][j])/erdt,4)
            if type == 2:
                option_table[i][j] = round(max(option_table[i][j],iopt*(share_table[i][j]-X)),4)
    return option_table

def main():
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
            X = float(input('Strike Price:= '))
            if X < 0:
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
    while True:
        try:
            type = int(input('type -> 1:= European Style | 2:= American Style: '))
            if type != 1 and type != 2:
                continue
            else:
                break
        except ValueError as e:
            print('1:= European Style | 2:= American Style')
            continue
    if type == 1:
        type_output = 'European'
    else:
        type_output = 'American'
    if iopt == 1:
        iopt_output = 'call'
    else:
        iopt_output = 'put'
    print('\nThe '+type_output+' style '+iopt_output+' option premium is: {0:>10s}'.format('£'+str(option_binomial(S,X,rate,div,t,T,sig,step,iopt,type)[0][0])))

#--------------------------------------------------
main()