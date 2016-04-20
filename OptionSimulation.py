import numpy as np

# parameter initialization
mu,sigma = 0.1, 0.1
r = 0.05
S0 = 100
K = 100

N = 52
T = 1.0
dt = T/N
n = 1000

# random number generation
NR = np.random.normal(0,1,[n,N])

#  Problem 1 (a) and (b)
ST = np.ones((1000,53)) # from S0 = 100 for all paths to S(52*dt)=S(1)
paths = range(1000)
timeStep = range(1,53,1)
Savg = np.zeros(1000)
payoff = np.zeros(1000)
STanti = np.ones((1000,53))
antipayoff = np.zeros(1000)
Savganti = np.zeros(1000)
Vpayoff = np.zeros(1000)

Gavg = np.ones(1000)
Gpayoff = np.zeros(1000)

for p in paths:
    ST[p,0] = S0
    STanti[p,0] = S0
    
    for t in timeStep:
        # generate ST(i*dt)
        ST[p,t] = ST[p,t-1]*np.exp((r-0.5*pow(sigma,2))*dt+sigma*np.sqrt(dt)*NR[p,t-1])
        # save the average price
        Savg[p] = (Savg[p]*(t-1) + ST[p,t])/t
        
        STanti[p,t] = STanti[p,t-1]*np.exp((r-0.5*pow(sigma,2))*dt+sigma*np.sqrt(dt)*(-1)*NR[p,t-1])
        Savganti[p] = (Savganti[p]*(t-1) + STanti[p,t])/t

        Gavg[p] = np.power(np.power(Gavg[p],t-1)*ST[p,t],1.0/t)

    # calculate discounted payoff
    payoff[p] = np.exp(-r*T)*max([(Savg[p] - K),0])
    antipayoff[p] = np.exp(-r*T)*max([(Savganti[p] - K),0])
    Gpayoff[p] = np.exp(-r*T)*max([(Gavg[p] - K),0])
    
    # antithetic estimate
    Vpayoff[p] = (payoff[p] + antipayoff[p])/2

# take an average over the paths and calculate standard error of estimates
price = np.average(payoff)
se = np.std(payoff)/np.sqrt(n)
print price
print se

antiprice = np.average(Vpayoff)
seanti = np.std(Vpayoff)/np.sqrt(n)
print antiprice
print seanti

# Problem 1 (c)
STav = np.average(ST[:,52])
y = payoff
x = ST[:,52]
tmp = np.corrcoef(x,y)
rho = tmp[0][1]
ahat = (-1)*rho*np.std(y)/np.std(x)

CVprice = price + ahat*(STav - S0*np.exp(r*T))
CVse = np.std(y)*np.sqrt(1-rho*rho)/np.sqrt(n)
print CVprice
print CVse

# Problem 1 (d)
from math import log, sqrt, exp
from scipy.stats import norm

def call_option_pricer(S0, K, T, r, sigma, d):
    
    d1 = (log(S0*exp(-d*T)/K) + (r + 0.5 * sigma *sigma) * T) / sigma / sqrt(T)
    d2 = d1 - sigma * sqrt(T)
    
    price = S0*exp(-d*T) * norm.cdf(d1) - K * exp(-r*T) * norm.cdf(d2)
    return price

d = 0
sigma_star = sigma/sqrt(3)
d_star = (r+d)/2 + sigma*sigma/12
CVex = call_option_pricer(S0,K,T,r,sigma_star,d_star)

CVav = np.average(Gpayoff)
y = payoff
x = Gpayoff
tmp = np.corrcoef(x,y)
rho = tmp[0][1]
ahat = (-1)*rho*np.std(y)/np.std(x)

CVprice = price + ahat*(CVav - CVex)
CVse = np.std(y)*np.sqrt(1-rho*rho)/np.sqrt(n)
print CVprice
print CVse

# Problem 1 (f)
d = 0
sigma_star = sigma*sqrt((N+1)*(2*N+1.0)/(6*N*N))
d_star = (r*(N-1.0)/N+d*(N+1.0)/N)/2 + sigma*sigma*(N*N-1.0)/(N*N)/12
CVex = call_option_pricer(S0,K,T,r,sigma_star,d_star)

CVav = np.average(Gpayoff)
y = payoff
x = Gpayoff
tmp = np.corrcoef(x,y)
rho = tmp[0][1]
ahat = (-1)*rho*np.std(y)/np.std(x)

CVprice = price + ahat*(CVav - CVex)
CVse = np.std(y)*np.sqrt(1-rho*rho)/np.sqrt(n)
print CVprice
print CVse

# Comment on the effectiveness of Variance Reduction Methods:
# 1. antithetic variables is pretty good, it reduces the variance significantly
# 2. the effectiveness of control variables depend largely on whether you can find an estimate highly correlated to what you want to estimate

# Problem 2 (a) and (b)

# parameter initialization
mu,sigma = 0.1, 0.2
r = 0.05
S0 = 100
K = [120.0,140.0,160.0]
T = 1.0
n = 10000

# random number generation
NR = np.random.normal(0,1,n)

ST = np.ones(n) 
paths = range(n)
payoff = np.zeros((n,3))
price = np.zeros(3);
se = np.zeros(3);
putpayoff = np.zeros((n,3))
putprice = np.zeros(3)
seput = np.zeros(3)

for i in range(3):
    for p in paths:
        ST[p] = S0*np.exp((r-0.5*sigma*sigma)*T+sigma*np.sqrt(T)*NR[p])
        payoff[p,i] = np.exp(-r*T)*max([(ST[p] - K[i]),0])
        putpayoff[p,i] = np.exp(-r*T)*(max([(K[i]-ST[p]),0]))
    price[i] = np.average(payoff[:,i])
    putprice[i] = np.average(putpayoff[:,i])+S0-K[i]*np.exp(-r*T)
    se[i] = np.std(payoff[:,i])/np.sqrt(n)
    seput[i] = np.std(putpayoff[:,i])/np.sqrt(n)

# Problem 2 (c)
CVav = np.average(ST)
CVex = S0*np.exp(r*T)
CVprice = np.zeros(3)
CVse = np.zeros(3)
for i in range(3):
    y = putpayoff[:,i]
    x = ST
    tmp = np.corrcoef(x,y)
    rho = tmp[0][1]
    ahat = (-1)*rho*np.std(y)/np.std(x)
    CVprice[i] = price[i] + ahat*(CVav - CVex)
    CVse[i] = np.std(y)*np.sqrt(1-rho*rho)/np.sqrt(n)
print CVprice
print CVse

ST = np.ones(n) 
payoff = np.zeros((n,3))
ISprice = np.zeros(3)
ISse = np.zeros(3)
for i in range(3):
    L = (log(K[i]/S0)-(r-0.5*sigma*sigma)*T)/sigma/np.sqrt(T)
    for p in paths:
        U = np.random.uniform();
        X = norm.ppf(U*(1-norm.cdf(L))+norm.cdf(L))
        ST[p] = S0*np.exp((r-0.5*sigma*sigma)*T+sigma*np.sqrt(T)*X)
        payoff[p,i] = np.exp(-r*T)*(ST[p]-K[i])*(1-norm.cdf(L))

    ISprice[i] = np.average(payoff[:,i])
    ISse[i] = np.std(payoff[:,i])/np.sqrt(n)

print ISprice
print ISse


# clearly we can see that Importance Sampling do a much better job than other estimates

# Problem 3
T = 0.25
m = 50
r = 0.05
sigma = 0.15
S0 = 95
n = 100000
dt = T/m

# random number generation
NR = np.random.normal(0,1,[n,m])

ST = np.ones((n,m+1))
paths = range(n)
timeStep = range(1,m+1,1)
payoff = np.zeros(n)
Smin = np.zeros(n)

H = [94,90,85,90]
K = [96,96,96,106]

price = np.zeros(4)
se = np.zeros(4)

for i in range(4):
    print i
    for p in paths:
        ST[p,0] = S0
        Smin[p] = 100000.0
    
        for t in timeStep:
            # generate ST(i*dt)
            ST[p,t] = ST[p,t-1]*np.exp((r-0.5*pow(sigma,2))*dt+sigma*np.sqrt(dt)*NR[p,t-1])
            Smin[p] = np.min([Smin[p],ST[p,1],ST[p,t]])
        
        if Smin[p]<H[i] and ST[p,t]>K[i]:
            payoff[p] = np.exp(-r*T)
        else:
            payoff[p] = 0

    price[i] = np.average(payoff)
    se[i] = np.std(payoff)/np.sqrt(n)

print price
print se
# Problem 3 (b)

# i
def digital_option_pricer(S0, K, T, r, sigma, d):
    if T==0:
        if (S0>K):
            return 1
        else:
            return 0
    d1 = (log(S0*exp(-d*T)/K) + (r + 0.5 * sigma *sigma) * T) / sigma / sqrt(T)
    d2 = d1 - sigma * sqrt(T)
    
    price = exp(-r*T) * norm.cdf(d2)
    return price

# ii

ST = np.ones((n,m+1))
paths = range(n)
timeStep = range(1,m+1,1)
payoff = np.zeros(n)
Smin = np.zeros(n)

H = [94.0,90.0,85.0,90.0]
K = [96.0,96.0,96.0,106.0]

priceCMC = np.zeros(4)
seCMC = np.zeros(4)

for i in range(4):
    print i
    for p in paths:
        ST[p,0] = S0
    
        for t in timeStep:
            # generate ST(i*dt)
            ST[p,t] = ST[p,t-1]*np.exp((r-0.5*pow(sigma,2))*dt+sigma*np.sqrt(dt)*NR[p,t-1])        
            
            if (ST[p,t]<H[i]):
                payoff[p] = exp(-r*t*dt)*digital_option_pricer(ST[p,t], K[i], T-t*dt, r, sigma, 0)
                break

    priceCMC[i] = np.average(payoff)
    seCMC[i] = np.std(payoff)/np.sqrt(n)

print priceCMC
print seCMC

vr = se**2/seCMC**2
# We should be more confident with the result given by Conditional Monte Carlo
# Although normal Monte Carlo seems to give a smaller standard error, it is very likely due to the fact that most paths aren't knocked-in

# Problem 4
T = 0.2
m = 50 # m = 25/50
r = 0.1
sigma = 0.3
S0 = 100
n = 100000
dt = T/m

# random number generation
NR = np.random.normal(0,1,[n,m])

ST = np.ones((n,m+1))
paths = range(n)
timeStep = range(1,m+1,1)
payoff = np.zeros(n)
Smin = np.zeros(n)

H = 95
K = 100

for p in paths:
        ST[p,0] = S0
        Smin[p] = 100000.0
    
        for t in timeStep:
            # generate ST(i*dt)
            ST[p,t] = ST[p,t-1]*np.exp((r-0.5*pow(sigma,2))*dt+sigma*np.sqrt(dt)*NR[p,t-1])
            Smin[p] = np.min([Smin[p],ST[p,1],ST[p,t]])
        
        if Smin[p]<H and ST[p,t]>K:
            payoff[p] = np.exp(-r*T)*(ST[p,t]-K)
        else:
            payoff[p] = 0

        price = np.average(payoff)
        se = np.std(payoff)/np.sqrt(n)

print price
print se

def down_and_in_call(S0, K, T, r, sigma, d,H):
    l = (r-d+sigma*sigma/2)/sigma**2
    y = log((H*1.0)**2/(S0*1.0*K))/(sigma*sqrt(T)) + l*sigma*sqrt(T)
   
    price = S0*exp(-d*T)*np.power(H*1.0/S0,2*l)*norm.cdf(y) - K*exp(-r*T)*np.power(H*1.0/S0,2*l-2) * norm.cdf(y-sigma*sqrt(T))
    return price

down_and_in_call(S0,K,T,r,sigma,0,H)

#spot = 2.45
#strike = 2.50
#maturity = 0.25
#r = 0.05
#vol = 0.25


#print 'Option Price : %.4f' % call_option_pricer(spot, strike, maturity, r, vol)

#import time


#portfolioSize = range(1, 1000, 50)
#timeSpent = []

#for size in portfolioSize:
#    now = time.time()
#    strikes = np.linspace(2.0,3.0,size)
#    for i in range(size):
#        res = call_option_pricer(spot, strikes[i], maturity, r, vol)
#        #print 'Option Price : %.4f' % res
#    timeSpent.append(time.time() - now)

#import matplotlib.pyplot as plt
#plt.plot(portfolioSize,timeSpent)
#plt.ylabel('T(Seconds)',fontsize = 15)
##plt.xlabel('Smarts')
#plt.title('Option Pricing Time',fontsize = 20)
##plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
#plt.axis([0, 1000, 0, 0.2])
#plt.grid(True)
#plt.show()
