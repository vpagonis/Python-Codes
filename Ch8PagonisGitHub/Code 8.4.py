# Plots of KV-LM equation for different model parameter
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
from scipy.special import wrightomega
def KVLM(R,Aopt,j):
    c=(n0/N)*(1-R)/R
    zLMOSL=(1/c)-np.log(c)+Aopt*(t**2.0)/(2*P)
    lam=wrightomega(zLMOSL)
    plt.plot(t,(N*R/((1-R)**2.0))*Aopt*t/(lam+lam**2),symbs[j],
    label=labls[j])
N,     R,  Aopt, n0 =\
1e10, 0.1, 0.3,  1e9
t = np.linspace(0, 100, 100)
P=max(t)
plt.subplot(1,2, 1)
symbs=['o-','+-','^-']
labls=['R=0.1','R=0.2','R=0.3']
KVLM(0.1,0.3,0)
KVLM(0.2,0.3,1)
KVLM(0.3,0.3,2)
plt.title('(a)')
plt.text(15,.3e9,'KV-LM equation')
plt.text(15,.15e9,'for LM-OSL')
plt.text(15,.01e9,'Variable R')
plt.xlabel('Time [s]')
plt.ylabel('LM-OSL [a.u.]')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.subplot(1,2, 2)
labls=[r'$\lambda$=0.3 s$^{-1}$',r'$\lambda$=0.4 s$^{-1}$',\
r'$\lambda$=0.5 s$^{-1}$']
KVLM(0.1,.3,0)
KVLM(0.1,.4,1)
KVLM(0.1,.5,2)
plt.title('(b)')
plt.text(7,.3e9,'KV-LM equation')
plt.text(7,.1e9,'Variable power')
plt.xlabel('Time [s]')
plt.ylabel('LM-OSL [a.u.]')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.tight_layout()
plt.show()