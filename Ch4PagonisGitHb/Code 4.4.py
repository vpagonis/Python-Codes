import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
# Simulate Loss of charge due to ground state tunneling (GST)
z=1.8                   # constant in GST model
s=3e15                  # frequency factor
years=1000              # number of years elapsed
elapsedt=3.154e7*years  # change years to seconds
times=np.arange(1,elapsedt,1e6)          # times in seconds
rhos=[1e-6,5e-6,1e-5]                    # rho-prime values
##### function to find n(t) ###
def n(rho):
    return 100*np.exp(-rho*(np.log(z*s*times)**3.0))
ns=[n(x) for x in rhos] 
plt.plot(times/3.154e7,ns[0],'-',c='black',\
linewidth=2,label=r'$\rho$'+"'"+'=10$^{-6}$')
plt.plot(times/3.154e7,ns[1],'-.',c='r',\
linewidth=2,label=r'$\rho$'+"'"+'=5x10$^{-6}$')
plt.plot(times/3.154e7,ns[2],'--',c='b',\
linewidth=2,label=r'$\rho$'+"'"+'=10$^{-5}$')
leg = plt.legend()
plt.ylim(0,120)
leg.get_frame().set_linewidth(0.0)
plt.xlabel('Elapsed time t [a]')
plt.ylabel('Remaining electrons [%]')
plt.title('Loss of charge in ground state tunneling')
plt.tight_layout()
plt.show() 