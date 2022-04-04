# Plot of the analytical equations for MOK
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
kB,       E,  s,    beta =\
8.617e-5, 1, 1e12,   1 
def nMOK(tmps):
    expint=kB*((273+tmps)**2.0)/E*\
	    np.exp(-E/(kB*(273+tmps)))*(1-2*kB*(273+tmps)/E)
    Ft=np.exp(g*Nd*s/beta*expint)
    return alpha*Nd/(Ft-alpha)
def TLMOK(tmps):
    expint=kB*((273+tmps)**2.0)/E*\
	    np.exp(-E/(kB*(273+tmps)))*(1-2*kB*(273+tmps)/E)
    Ft=np.exp(g*Nd*s/beta*expint)
    return g*(Nd**2.0)*alpha*s*np.exp(-E/(kB*(273+tmps)))*\
    Ft/((Ft-alpha)**2.0)
tims= range(25,160)
temps=[beta*tim for tim in tims]
n10, N1, Nd= 1e10 ,1e10, 1e12
alpha, g =n10/(n10+Nd), 1/(N1+Nd)
plt.subplot(1,2, 1)
plt.plot(temps,[nMOK(x) for x in temps],'o',label='n(t)')
plt.plot(temps,[30*TLMOK(x) for x in temps],'+',
label='TL x30')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.title('(a)')
plt.ylabel('n(t) and TL (MOK)')
plt.xlabel(r'Temperature T [$^{o}$C]')
plt.text(25,.8e10,'MOK model')
plt.subplot(1,2, 2)
plt.plot(temps,[TLMOK(x) for x in temps],'+-',
label=r'$\alpha$=0.01')
n10, N1, Nd= 1e10 ,1e10, 1e10
alpha, g =n10/(n10+Nd), 1/(N1+Nd)
plt.plot(temps,[TLMOK(x) for x in temps],'^-',
label=r'$\alpha$=0.5')
n10, N1, Nd= 1e10 ,1e10, 1e8
alpha, g =n10/(n10+Nd), 1/(N1+Nd)
plt.plot(temps,[TLMOK(x) for x in temps],'o-',
label=r'$\alpha$=0.99')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.title('(b)')
plt.xlim(40,160)
plt.ylabel('TL [a.u.]')
plt.xlabel(r'Temperature T [$^{o}$C]')
plt.tight_layout()
plt.show()