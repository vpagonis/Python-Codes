#Example of summing partial TL curves simplest
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
##### function to find distribution of distances ###
def partialTL(tmps,rprime):
    seff=s*np.exp(-(rho**(-1/3))*rprime)
    TLpartial=3*(rprime**2.0)*np.exp(-(rprime**3.0))*\
	    (seff/hr)*np.exp(-E/(kB*(273+tmps)))*\
    np.exp(-seff*rprime*kB*((273+tmps)**2.0)/(hr*E)*\
    np.exp(-E/(kB*(273+tmps)))*(1-2*kB*(273+tmps)/E))
    return TLpartial*dr
kB,         s,        E ,    hr=\
8.617e-5,  2e15,  1.3 ,   1
t = np.linspace(0, 500, 500)
temps=hr*t
rho=1e-3                       # rho-prime value
dr=0.1
rprimes=np.arange(0,2.2,dr)    # rprime=0-2.2
plt.subplot(1,2, 1)
for i in range(len(rprimes)):
    plt.plot(temps,[partialTL(x,rprimes[i]) for x in temps])
plt.text(250,.0035,"Partial TL for ")
plt.text(250,.0032,"r'=0-2.2")
plt.ylabel('Partial TL signal [a.u.]')
plt.xlabel(r'Temperature [$^{o}$C]')
plt.title('(a)')
plt.subplot(1,2,2)
u=np.array([[partialTL(x,rprimes[i]) for x in temps]\
for i in range(len(rprimes))])
plt.plot(temps,sum(u),c='b')
plt.text(300,.01, 'Sum of')
plt.text(300,.009, 'partial')
plt.text(300,.008, 'TL curves')
plt.ylabel('Remnant TL signal [a.u.]')
plt.xlabel(r'Temperature [$^{o}$C]')
plt.title('(b)')
plt.tight_layout()
plt.show()