# Varying the parameter R in the OTOR model
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
    sol = solve_ivp(lambda t, y: f(t, y), [50,200],yinit)
    nvals,ncvals,mvals=np.squeeze(sol.y)
    temps=hr*np.squeeze(sol.t)
    return temps, nvals,ncvals,mvals
#######
temps, nvals, ncvals, mvals=TL(1e-11,1,1e12)
plt.plot(temps, Am*ncvals*mvals,'-',linewidth=3,
label='R=0.001')
temps, nvals, ncvals, mvals=TL(1e-8,1,1e12)
plt.plot(temps, Am*ncvals*mvals,'--',linewidth=3,
label='R=1')
temps, nvals, ncvals, mvals=TL(3e-8,1,1e12)
plt.plot(temps, Am*ncvals*mvals,'-.',linewidth=3,
 label='R=3')
plt.ylim(0,4e8)
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('TL signal [a.u.]')
plt.xlabel(r'Temperature T [$^{o}$C]')
plt.title(r'OTOR model, varying the retrapping ratio R')
plt.text(170,2.5e8,r'R=A$_{n}$/A$_{m}$')
plt.tight_layout()
plt.show()