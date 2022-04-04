#plot KP-LM equation for localized processes
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
# pLM in this code represents the effective frequency s'  
def KPLM(rho,pLM):    
    LM=np.exp (-rho*(np.log(1 + pLM*(t**2/(2*P))))\
	    ** 3.0)*(np.log(1+pLM*(t**2/(2*P)))**2.0)/(1+(t**2/(2*P))*pLM)  
    plt.plot(t,LM/max(LM),symbs[j-1], linewidth=2,
    label=labls[j-1])
t = np.linspace(1, 100, 100)
P=max(t)
pLM=3
rhos=[1e-3,5e-3,1e-2]
labls=[r'$\rho$'+"'="+str(x) for x in rhos]
symbs=['+-','^-','o-']
plt.subplot(1,2,1) 
for j in range(1,4): 
    KPLM(rhos[j-1],pLM)
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.xlabel('Time [s]')
plt.text(50,.75,'KP-LM equation')
plt.ylabel('LM-IRSL [a.u]')
plt.title('(a)')
plt.subplot(1,2,2) 
pLM=[3,5,7]
labls=["s'="+str(x)+r' s$^-1$' for x in pLM]
for j in range(1,4): 
    KPLM(0.01,pLM[j-1])
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.xlabel('Time [s]')
plt.ylabel('LM-IRSL [a.u]')
plt.title('(b)') 
plt.tight_layout()
plt.show()