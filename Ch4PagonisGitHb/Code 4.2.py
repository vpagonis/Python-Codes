import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
s=3e15                 # frequency factor
rho=1e-6               # rho-prime values 0.005-0.02
rc=0.0                 # for freshly irradiated samples, rc=0
times=(0,3.154e9,3.154e11,3.154e13)       # times in seconds
rprimes=np.arange(0,2.2,0.002)            # rprime=0-2.2
##### function to find distribution of distances ###
def findDistr(tim):
    return 3*(rprimes**2.0)*np.exp(-(rprimes**3.0))*\
    np.exp(-np.exp(-(rho**(-1/3))*rprimes)*s*tim)
distribs=[findDistr(x) for x in times]
plt.plot(rprimes,distribs[0],'-',c='black',linewidth=4,\
label=r't=0')
plt.plot(rprimes,distribs[1],'+-',c='r',label=r'10$^{2}$ a')
plt.plot(rprimes,distribs[2],'--',c='b',label=r'10$^{4}$ a')
plt.plot(rprimes,distribs[3],label=r'10$^{6}$ a')
plt.ylim(0,1.3)
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('Nearest neighbor Distribution g(r$^{\'}$)')
plt.xlabel('Dimensionless distance r$^{\'}$')
plt.title('Time evolution of nearest neighbor distribution')
plt.tight_layout()
plt.show()