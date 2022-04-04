#plot KV-ITL equation for delocalized processes
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
from scipy.special import wrightomega
def KVITL(R,pITL):
    c=(n0/N)*(1-R)/R
    zITL=(1/c)-np.log(c)+(pITL*n0/(c*N*R))*t
    lam=wrightomega(zITL) 
    ITL=(N*R/((1-R)**2.0))*pITL/(lam+lam**2)   
    plt.plot(t,ITL/max(ITL),symbs[j-1],\
    linewidth=2,label=labls[j-1])
N,     s,   E,  n0 ,  kB=\
1e10, 1e12, 1,  1e9,  8.617e-5
t = np.linspace(0, 100, 100)
Tiso=110+273.15
pITL=s*np.exp(-E/(kB*Tiso))
Rs=[0.01,.1,.95]
symbs=['+','^','o']
labls=['R='+str(x) for x in Rs]
plt.subplot(1,2,1) 
for j in range(1,4): 
    KVITL(Rs[j-1],pITL)
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.xlabel('Time [s]')
plt.text(45,.7,'KV-ITL equation')
plt.ylabel('ITL [a.u]')
plt.title('(a)') 
plt.subplot(1,2,2) 
Tiso=[110,120,130] 
pITL=[s*np.exp(-E/(kB*(x+273.15))) for x in Tiso]
labls=['T='+str(x)+r'$^o$C' for x in Tiso]
for j in range(1,4): 
    KVITL(0.1,pITL[j-1])
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.xlabel('Time [s]')
plt.ylabel('ITL [a.u]')
plt.title('(b)') 
plt.tight_layout()
plt.show()