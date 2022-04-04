#Simulation of remnant TL after fading due to GST in nature
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
##### function to find distribution of distances in GST model
### after fading in nature for time tim
def findDistr(tim):
    return 3*(rprimes**2.0)*np.exp(-(rprimes**3.0))*\
    np.exp(-np.exp(-(rho**(-1/3))*rprimes)*stun*tim)
# find partial TL using the EST model
def partialTL(tmps,rprime,tim):
    seff=sEST*np.exp(-(rhoEST**(-1/3))*rprime)
    distr=3*(rprime**2.0)*np.exp(-(rprime**3.0))*\
    np.exp(-np.exp(-(rho**(-1/3))*rprime)*stun*tim)
    TLpartial=distr*\
	    (seff/hr)*np.exp(-E/(kB*(273+tmps)))*\
    np.exp(-seff*kB*((273+tmps)**2.0)/(hr*E)*\
    np.exp(-E/(kB*(273+tmps)))*(1-2*kB*(273+tmps)/E))
    return TLpartial*dr
# Here we use s and rho values for GST model
stun=3e15                  			# frequency factor
rho=3e-6               				# rho-prime value
timesNature=(3600*24*365,100*3600*24*365,1e4*3600*24*365)
# fading times in nature, in seconds
dr=0.02
rprimes=np.arange(0,2.2,dr)    	 # rprime=0-2.2
styls=['solid','dotted','dashed']
labls=[r't=1 a',r'10$^{2}$ a',r'10$^{4}$ a']
plt.subplot(1,2, 1)
for j in range(3):
    plt.plot(rprimes,findDistr(timesNature[j]),linestyle=\
    styls[j],label=labls[j])
plt.ylim(0,1.3)
plt.xlim(0,2.3)
plt.ylabel('Trapped electron distribution n(r$^{\'}$,t)')
plt.xlabel('Dimensionless distance r$^{\'}$')
plt.title('(a)')
leg = plt.legend(loc='upper right')
leg.get_frame().set_linewidth(0.0)
## use different rhoEST and sEST for the EST model
kB,         sEST,        E ,    rhoEST,   hr=\
8.617e-5,  3.5e12,  1.45 ,    5e-3,   1
t = np.linspace(0, 500, 500)
temps=hr*t
plt.subplot(1,2, 2) 
plt.title('(b)')
for j in range(3):
    u=np.array([[partialTL(x,rprimes[i],timesNature[j]) for\
    x in temps] for i in range(len(rprimes))])
    plt.plot(temps,sum(u),linestyle=styls[j],label=labls[j])
plt.ylabel('Remnant TL signal [a.u.]')
plt.xlabel(r'Temperature [$^{o}$C]')
plt.xlim(150,500)
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.tight_layout()
plt.show()