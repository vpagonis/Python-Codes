# Numerical solution of the GOT equation for TL
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
# Numerical integration of GOT model.
n0,    N,   R,   s,   E ,hr , kB=\
1e10, 1e10, 1,  1e12, 1 , 1,  8.617e-5
t = np.linspace(0, 180, 180)
def deriv(y, t):
    n = y
    dndt = - (n**2.0)*s*np.exp(-E/(kB*(273+hr*t)))/((R*(N-n)+n))
    return dndt
y0 = n0
ret = odeint(deriv, y0, t)
n= ret.flatten()
# Plot the data 
plt.plot(hr*t, n, 'b+',  linewidth=3, label='n(t)')
plt.ylabel('n(t)  and  TL')
plt.ylim(0,1.2e10)
plt.xlim(0,180)
plt.title('Numerical integration of GOT model')
plt.xlabel(r'Temperature T [$^{o}$C]')
TL=  (n**2.0)*s*np.exp(-E/(kB*(273+hr*t)))/((R*(N-n)+n))
plt.plot(hr*t,30*TL,'r-', linewidth=3,label='TL x30')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.tight_layout()
plt.show()