# Numerically solve the ODE for R-W: first order kinetics
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
# Function for analytical solution
def nfirst(n0,tmps):
    return n0*np.exp(-s*kB*((273+tmps)**2.0)/(hr*E)*\
    np.exp(-E/(kB*(273+tmps)))*(1-2*kB*(273+tmps)/E))
def TLfirst(n0,tmps):
    return n0*s*np.exp(-E/(kB*(273+tmps)))*\
    np.exp(-s*kB*((273+tmps)**2.0)/(hr*E)*\
    np.exp(-E/(kB*(273+tmps)))*(1-2*kB*(273+tmps)/E))
# Numerical parameters for R-W model.
kB,        N,     s,    E , n0 ,  hr=\
8.617e-5, 1e10,  1e12,  1 , 1e10, 1
# A grid of time points and temperatures
t = np.linspace(0, 180, 60)
temps=hr*t
# The R-W model differential equation.
def deriv(y, t):
    n = y
    dndt = - n*s*np.exp(-E/(kB*(273+hr*t)))
    return dndt
y0 = n0
ret = odeint(deriv, y0, t)
n= ret.flatten()
# Plot the data 
plt.subplot(1,2, 1)
plt.plot(temps, n, 'o', label='Numerical')
plt.plot(temps,[nfirst(n0,x) for x in temps],'-',\
label='Analytical')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('Filled traps [cm$^{-3}$]')
plt.ylim(0,1.2e10)
plt.xlim(0,180)
plt.title('(a)')
plt.text(20,.6e10,'R-W model')
plt.text(20,.5e10,'numerical')
plt.text(20,.4e10,'solution n(t)')
plt.xlabel(r'Temperature T [$^{o}$C]')
plt.subplot(1,2, 2)
TL= n*s*np.exp(-E/(kB*(273+hr*t)))
plt.plot(temps,TL,'o',label='Numerical')
plt.plot(temps,[TLfirst(n0,x) for x in temps],'-',
label='Analytical')
plt.ylim(0,4e8)
plt.text(50,1.5,'TL')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('TL signal [a.u.]')
plt.xlabel(r'Temperature T [$^{o}$C]')
plt.title('(b)')
plt.tight_layout()
plt.show()