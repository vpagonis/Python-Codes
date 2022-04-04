# GOT MODEL- code for IRRADIATION Lambert solution
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
from scipy.special import lambertw
# Numerical parameters for GOT irradiation model.
N,     R, X,    n0 =\
1E10, .1, 1e10, 0
t = np.linspace(0, 100, 100)
# The GOT model irradiation differential equation.
def deriv(y, t):
    n = y
    dndt=   (N-n) * X*R / ((N - n) * R + n )
    return dndt
y0 = n0
ret = odeint(deriv, y0, t)
n= ret.flatten()
plt.plot(t, n, 'r^',  lw=2, label='GOT model')
ret = odeint(deriv, y0, t)
#analytical n(t)
Dc=N/R
u=N*(1+np.real(lambertw((R-1)*np.exp((R-1)-\
(t*X/np.abs(Dc))))/(1-R)))
plt.plot(t, u,label='PKC equation')
leg = plt.legend(loc='right')
leg.get_frame().set_linewidth(0.0)
plt.ylabel('Trap filling ratio, n/N')
plt.ylim(0,1.2e10)
plt.xlim(0,100)
plt.title('Dose response in GOT model, with PKC eqt')
plt.text(20,.6e10,'GOT model')
plt.text(20,.5e10,'Analytical solution')
plt.text(20,.4e10,'with Lambert')
plt.xlabel('Time [s]')
plt.tight_layout()
plt.show()