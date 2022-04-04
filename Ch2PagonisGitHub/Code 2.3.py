# Varying the parameters E,s in the OTOR model
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
# Numerical parameters for OTOR.
N,      Am,     hr ,  kB= \
1e10,  1e-8   , 1,  8.617e-5
yinit = [1e10,0,1e10]
def TL(An,E,s):
    def f(t, y):
        dydt =[-y[0]*s*np.exp(-E/(kB*(273+hr*t)))+y[1]*An*\
        (N-y[0]), y[0]*s*np.exp(-E/(kB*(273+hr*t)))- y[1]*\
        An*(N-y[0])-y[2]*Am*y[1], -y[2]*Am*y[1]]
        return dydt
    sol = solve_ivp(lambda t, y: f(t, y), [20,220],yinit)
    nvals,ncvals,mvals=np.squeeze(sol.y)
    temps=hr*np.squeeze(sol.t)
    return temps, nvals,ncvals,mvals
plt.subplot(1,2,1)
labels=['E=1.0 eV','E=1.05 eV','E=1.1 eV']
lins=['solid','dashed','dotted']
E=[1.,1.05,1.1]
for j in range (0,3,1): 
     temps, nvals, ncvals, mvals=TL(1e-10,E[j],1e12)
     plt.plot(temps, Am*ncvals*mvals,
     linestyle=lins[j], linewidth=3,label=labels[j])
plt.ylim(0,4.2e8)
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('TL signal [a.u.]')
plt.xlabel(r'Temperature T [$^{o}$C]')
plt.title(r'(a)')
plt.subplot(1,2,2)
s=[1e11,1e12,5e12]
labels=[r's=10$^{11}$ s$^{-1}$',r's=10$^{12}$ s$^{-1}$',
r's=5x10$^{12}$ s$^{-1}$']
for j in range (0,3,1): 
     temps, nvals, ncvals, mvals=TL(1e-10,1,s[j])
     plt.plot(temps, Am*ncvals*mvals,
     linestyle=lins[j],linewidth=3, label=labels[j])
plt.ylim(0,5.2e8)
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('TL signal [a.u.]')
plt.xlabel(r'Temperature T [$^{o}$C]')
plt.title(r'(b)')
plt.tight_layout()
plt.show()