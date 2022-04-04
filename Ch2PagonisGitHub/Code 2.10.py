# Plot the same TL glow curve for three different heating rates
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
kB,       E,  s,     n0 ,     N=\
8.617e-5, 1, 1e12,   1e10 ,   1e10
def TLfirst(beta,tmps):
    return n0*(s/beta)*np.exp(-E/(kB*(273+tmps)))*\
    np.exp(-s*kB*((273+tmps)**2.0)/(beta*E)*\
    np.exp(-E/(kB*(273+tmps)))*(1-2*kB*(273+tmps)/E))
temps=range(50,160,1)
symb=('o-','^-','+-')
for j in range(1,4,1):
    plt.plot(temps,[TLfirst(j,x) for x in temps],
    symb[j-1])
plt.title('Variable heating rates')
plt.xlim(50,160)
plt.ylim(0,4e8)
plt.text(105,3.3e8,'1 K/s')
plt.text(117,3.1e8,'2 K/s')
plt.text(130,2.8e8,'3 K/s')
plt.ylabel(r'TL/$\beta$ [a.u.]')
plt.xlabel(r'Temperature T [$^{o}$C]')
plt.show()