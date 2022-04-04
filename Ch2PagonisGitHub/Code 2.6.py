# Simulate first-order glow peaks with various 
# initial electron trap concentrations (n0).
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
kB,       E,  s,    beta =\
8.617e-5, 1, 1e12, 1
def TLfirst(n0,tmps):
    return n0*s*np.exp(-E/(kB*(273+tmps)))*\
    np.exp(-s*kB*((273+tmps)**2.0)/(beta*E)*\
    np.exp(-E/(kB*(273+tmps)))*(1-2*kB*(273+tmps)/E))
tims= range(50,int(170/beta)-1,1)
temps=[beta*tim for tim in tims]
plt.plot(temps,[TLfirst(1e10,x) for x in temps],'+-',c='r',
label=r'n$_{0}$=1e10 cm$^{-3}$')
plt.plot(temps,[TLfirst(2e10,x) for x in temps],'o-',c='g',
label=r'n$_{0}$=2e10 cm$^{-3}$')
plt.plot(temps,[TLfirst(3e10,x) for x in temps],'^-',c='b',
label=r'n$_{0}$=3e10 cm$^{-3}$')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.title('Variable Initial occupancy n0')
plt.text(60,8e8,'First order')
plt.text(60,7.2e8,'Kinetics')
plt.ylabel(r'TL [a.u.]')
plt.xlabel(r'Temperature T [$^{o}$C]')
plt.show()
