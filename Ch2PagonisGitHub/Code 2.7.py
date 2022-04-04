# Simulate second-order glow peaks with various 
# initial electron trap concentrations (n0).
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
kB,       E,  s,    beta , N=\
8.617e-5, 1, 1e12,   1 ,   1e10
def TLsec(n0,tmps):
    expint=kB*((273+tmps)**2.0)/E*\
	    np.exp(-E/(kB*(273+tmps)))*(1-2*kB*(273+tmps)/E)
    return (n0**2.0)*(s/N)*np.exp(-E/(kB*(273+tmps)))*\
    ((1+(n0*s/(beta*N))*expint)**(-2.00) )
tims= range(50,int(220/beta)-1,1)
temps=[beta*tim for tim in tims]
plt.plot(temps,[TLsec(2e9,x) for x in temps] ,'+-',
c='r',label=r'n$_{0}$/N=0.2')
plt.plot(temps,[TLsec(4e9,x) for x in temps] ,'o-',
c='g',label=r'n$_{0}$/N=0.4')
plt.plot(temps,[TLsec(6e9,x) for x in temps] ,'^-',
c='b',label=r'n$_{0}$/N=0.6')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.title('Variable Initial occupancy n0')
plt.text(50,1e8,'Second order')
plt.text(50,.9e8,'Kinetics')
plt.ylabel(r'TL [a.u.]')
plt.xlabel(r'Temperature T [$^{o}$C]')
plt.show()