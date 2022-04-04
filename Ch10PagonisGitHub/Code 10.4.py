#plot KP-CW equation for localized processes
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
# pCW in this code represents the effective frequency s'  
def KPCW(rho,pCW):    
    CW=np.exp (-rho*(np.log(1 + pCW*t))\
	    ** 3.0)*(np.log(1+pCW*t)**2.0)/(1+t*pCW)  
    plt.plot(t,CW/max(CW),symbs[j-1], linewidth=2,
    label=labls[j-1])
t = np.linspace(1, 100, 100)
pCW=3
rhos=[1e-3,5e-3,1e-2]
labls=[r'$\rho$'+"'="+str(x) for x in rhos]
symbs=['+-','^-','o-']
plt.subplot(1,2,1) 
for j in range(1,4): 
    KPCW(rhos[j-1],pCW)
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.xlabel('Time [s]')
plt.text(45,.6,'KP-CW equation')
plt.ylabel('CW [a.u]')
plt.title('(a)')
plt.subplot(1,2,2) 
pCW=[3,5,7]
labls=["s'="+str(x)+r' s$^-1$' for x in pCW]
for j in range(1,4): 
    KPCW(0.01,pCW[j-1])
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.xlabel('Time [s]')
plt.ylabel('CW [a.u]')
plt.title('(b)') 
plt.tight_layout()
plt.show()