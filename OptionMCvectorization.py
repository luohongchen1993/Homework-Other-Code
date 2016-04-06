# vectorization trial No. 1
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
NR = np.random.normal(0,1,[n,N+1])

#  Problem 1 (a) and (b)
ST = np.ones([1000,53]) # from S0 = 100 for all paths to S(52*dt)=S(1)
paths = range(1000)
timeStep = range(1,53,1)
Savg = np.zeros([1000,1])

ST[:,0] = S0
for t in timeStep:
    # generate ST(i*dt)
    ST[:,t] = ST[:,t-1]*np.exp((r-0.5*pow(sigma,2))*dt+sigma*np.sqrt(dt)*NR[:,t-1])

Savg = np.average(ST,axis = 1)
payoff = np.exp(-r*T)*np.maximum(Savg-K,0)
print np.average(payoff)
print np.std(payoff)/np.sqrt(n)
