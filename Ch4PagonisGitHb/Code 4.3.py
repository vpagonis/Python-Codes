#Example of summing exponentials to obtain n(t)in GST
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
### function to find n(r',t) at time t and at distance r'###
def nrprimet(tim,rprime):
    seff=s*np.exp(-(rho**(-1/3))*rprime)
    nt=3*(rprime**2.0)*np.exp(-(rprime**3.0))*\
    np.exp(-seff*tim)
    return nt*dr
s=3e15                          # GST 
dr=0.01                         # interval in r'
rprimes=np.arange(0,2.2,dr)     # rprime=0-2.2
y=24*3600*365                   # year in s
times=10e3*y*np.linspace(0.01,2,30)
rho=1e-6                        # acceptor density
u=np.array([[nrprimet(x,rprimes[i]) for x in times]\
for i in range(len(rprimes))])
timesy=[i /y for i in times]
plt.plot(timesy,sum(u),'o',c='b',label=r'$\rho$'+"'"+'=\
1x10$^{-6}$')
ts=y*np.arange(0,2e4,10)
plt.plot(ts/y,np.exp(-rho*np.log(1.8*s*ts)**3.0))
plt.ylim(0.5,1)
plt.ylabel('Remaining charge  n(t)/no')
plt.xlabel('Elapsed time t [a]')
plt.title('Loss of charge in GST model')
rho=1.5e-6                        # acceptor density
u=np.array([[nrprimet(x,rprimes[i]) for x in times]\
for i in range(len(rprimes))])
plt.plot(timesy,sum(u),'^',c='r',label=\
r'$\rho$'+"'"+'=1.5x10$^{-6}$')
plt.plot(ts/y,np.exp(-rho*np.log(1.8*s*ts)**3.0))
leg = plt.legend()
plt.text(.95e4,.95,'Acceptor density')
leg.get_frame().set_linewidth(0.0)
plt.tight_layout()
plt.show()