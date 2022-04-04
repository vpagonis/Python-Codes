#  Plot the solution of GOT, using the Wright omega function 
from scipy.special import lambertw
from scipy.special import wrightomega
import numpy as np
import matplotlib.pyplot as plt
t = np.arange(0.0, 200.0, 1.0)
ys = lambertw(t)
no,      N,      R,    s,   E, kB,      beta= \
1.0E10, 1.0E10, .99, 1e12, 1, 8.617e-5, 1
c=(no/N)*(1-R)/R
expint=kB*((273+beta*t)**2.0)/(beta*E)*\
    np.exp(-E/(kB*(273+t*beta)))*(1-2*kB*(273+beta*t)/E)
zTL=(1/c)-np.log(c)+(s*expint/(beta*(1-R)))
lam=np.real(lambertw(np.exp(zTL)))
lam2=wrightomega(zTL)
plt.subplot(1,2, 1)
plt.plot(t,(N*R/((1-R)))*(1/lam),'r-',label='n(t)')
plt.plot(t,30*(N*R/((1-R)**2.0))*s*np.exp(-E/(kB*(273+t*beta)\
))/(lam+lam**2),'b--',label='TL x30')
plt.title('(a)')
plt.ylabel('TL signal [a.u.]')
plt.xlabel(r'Temperature T [$^{o}$C]')
plt.text(0,.7e10,'KV-TL eqt.')
plt.text(0,.6e10,'lambertw')
plt.text(0,.5e10,'overflow')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.subplot(1,2, 2)
plt.ylabel('n(t)   and  TL')
plt.text(0,.7e10,'KV-TL eqt.')
plt.text(0,.6e10,'with')
plt.text(0,.5e10,'wrightomega')
plt.xlabel(r'Temperature T [$^{o}$C]')
plt.plot(t,(N*R/((1-R)))*(1/lam2),'r-',label='n(t)')
plt.plot(t,30*(N*R/((1-R)**2.0))*s*np.exp(-E/(kB*(273+t*beta)\
))/(lam2+lam2**2),'b--',label='TL x30')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('TL signal [a.u.]')
plt.xlabel(r'Temperature T [$^{o}$C]')
plt.title('(b)')
plt.tight_layout()
plt.show()