#Example of summing exponentials to obtain n(t)in IGST model
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
### function to find n(r',t) at time t and at distance r'###
def nrprimet(tim,rprime,rho):
    seff=s*np.exp(-(rho**(-1/3))*rprime)
    expsoln=3*(rprime**2.0)*np.exp(-(rprime**3.0))*\
    X/(Do*seff+X)*(1-np.exp(-(X/Do+seff)*tim))
    return expsoln*dr
s=2e15                   		# frequency factor
Do=1600                      	# D0 in Gy
X=2.85/(1000*365*3600*24)    	# natural dose rate 2.85 Gy/ka
dr=0.1                       	# interval in r'
rprimes=np.arange(0,2.2,dr)     # rprime=0-2.2
y=24*3600*365                   # year in s
times=1e6*y*np.linspace(0,2,20)
timesy=[i /y for i in times]
rhos=[1e-6,2e-6,3e-6]          	# rho-prime values 0.005-0.02
styls=['o','^','+']
for j in range(3):
    u=np.array([[nrprimet(x,rprimes[i],rhos[j]) for x in times]\
    for i in range(len(rprimes))])
    plt.plot(timesy,sum(u),styls[j],label=r'$\rho$='+str(rhos[j]))
    ts=y*np.arange(0,2.1e6,1e5)
    plt.plot(ts/y,(1-np.exp(-(X/Do)*ts))*np.exp(-rhos[j]*\
    np.log(Do*s/X)**3.0))
plt.ylim(0,1.2)
plt.ylabel('Percent of filled traps,  n(t)/N')
plt.xlabel('Elapsed time [a]')
plt.title('IGST model numerical integration')
leg = plt.legend()
plt.text(1.e6,1.12,'Acceptor density')
leg.get_frame().set_linewidth(0.0)
plt.tight_layout()
plt.show()