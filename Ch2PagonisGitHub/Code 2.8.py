# Compare 1st and 2nd order TL with the same parameters
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
kB,       E,  s,    beta , N=\
8.617e-5, 1, 1e12,   1 ,   1e10
def TLfirst(n0,tmps):
    return n0*s*np.exp(-E/(kB*(273+tmps)))*\
    np.exp(-s*kB*((273+tmps)**2.0)/(beta*E)*\
    np.exp(-E/(kB*(273+tmps)))*(1-2*kB*(273+tmps)/E))
def TLsec(n0,tmps):
    expint=kB*((273+tmps)**2.0)/E*\
	    np.exp(-E/(kB*(273+tmps)))*(1-2*kB*(273+tmps)/E)
    return (n0**2.0)*(s/N)*np.exp(-E/(kB*(273+tmps)))*\
    ((1+(n0*s/(beta*N))*expint)**(-2.00) )
tims= range(25,int(220/beta)-1,1)
temps=[beta*tim for tim in tims]
plt.plot(temps,[TLfirst(1e10,x) for x in temps],'+-',
c='r',label=r'b=1')
plt.plot(temps,[TLsec(1e10,x) for x in temps] ,'o-',c='b',
label=r'b=2')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.title('Compare b=1 and b=2')
plt.text(40,2.7e8,'Compare Kinetics')
plt.text(40,2.5e8,'First order vs')
plt.text(40,2.3e8,'second order')
plt.ylabel(r'TL [a.u.]')
plt.xlabel(r'Temperature T [$^{o}$C]')
plt.show()