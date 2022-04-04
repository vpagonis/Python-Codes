#  Plot the analytical solution of GOT, using Lambert W-function 
from scipy.special import lambertw
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
t = np.arange(0.0, 150.0, 1.0)
ys = lambertw(t)
plt.subplot(1,2, 1)
plt.plot(t, ys,label='Lambert W(x)')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('W(x)')
plt.xlabel('x')
plt.title('(a)')
plt.subplot(1,2, 2)
t = np.arange(0.0, 150.0, 3)
no,      N,      R,    s,   E, kB,      beta= \
1.0E10, 1.0E10, 1e-4, 1e12, 1, 8.617e-5, 1
c=(no/N)*(1-R)/R
expint=kB*((273+beta*t)**2.0)/(beta*E)*\
    np.exp(-E/(kB*(273+t*beta)))*(1-2*kB*(273+beta*t)/E)
zTL=(1/c)-np.log(c)+(s*expint/(beta*(1-R)))
lam=np.real(lambertw(np.exp(zTL)))
plt.plot(t,(N*R/((1-R)))*(1/lam),'r-',label='n(t)')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('n(t)   and  TL')
plt.ylim(0,1.2e10)
plt.xlim(0,150)
plt.title('(b)')
plt.text(20,.6e10,'GOT model')
plt.text(20,.5e10,'KV-TL equation')
plt.xlabel(r'Temperature T [$^{o}$C]')
plt.plot(t,30*(N*R/((1-R)**2.0))*s*np.exp(-E/(kB*(273+t*beta)\
))/(lam+lam**2),'b--',label='TL x30')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('TL signal [a.u.]')
plt.xlabel(r'Temperature T [$^{o}$C]')
plt.title('(b)')
plt.tight_layout()
plt.show()