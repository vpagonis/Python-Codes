import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
from scipy.special import wrightomega
def KVCW(R,lamda,j):
    c=(n0/N)*(1-R)/R
    zCWOSL=(1/c)-np.log(c)+(lamda*n0/(c*N*R))*t
    lam=wrightomega(zCWOSL)
    plt.plot(t,(N*R/((1-R)**2.0))*(1/(lam+lam**2.0)),symbs[j],
    label=labls[j])
N,     n0 = 1e10,  1e10
t = np.linspace(0, 100, 100)
plt.subplot(1,2, 1)
symbs=['o-','+-','^-']
labls=['R=0.1','R=0.7','R=0.9']
KVCW(0.1,0.1,0)
KVCW(0.7,0.1,1)
KVCW(0.9,0.1,2)
plt.title('(a)')
plt.text(40,.6e10,'KV-CW equation')
plt.text(40,.5e10,'for CW-OSL')
plt.text(40,.4e10,'Variable R')
plt.xlabel('Time [s]')
plt.ylabel('CW-OSL [a.u.]')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.subplot(1,2, 2)
labls=[r'$\lambda$=0.1 s$^{-1}$',r'$\lambda$=0.2 s$^{-1}$',\
r'$\lambda$=0.3 s$^{-1}$']
KVCW(0.1,.1,0)
KVCW(0.1,.2,1)
KVCW(0.1,.3,2)
plt.title('(b)')
plt.text(40,.6e10,'KV-CW equation')
plt.text(40,.5e10,'Variable power')
plt.xlabel('Time [s]')
plt.ylabel('CW-OSL [a.u.]')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.tight_layout()
plt.show()