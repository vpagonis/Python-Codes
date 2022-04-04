#plot KP-ITL equation for localized processes
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")   
def KPITL(rho,pITL):    
    ITL=np.exp (-rho*(np.log(1 + z*pITL*t))\
	    ** 3.0)*(np.log(1+z*pITL*t)**2.0)/(1+t*pITL)  
    plt.plot(t,ITL/max(ITL),symbs[j-1], linewidth=2,
    label=labls[j-1])
s,   E,   kB,       z=\
1e12,  1,   8.617e-5, 1.8
t = np.linspace(1, 100, 100)
Tiso=170+273.15
pITL=s*np.exp(-E/(kB*Tiso))
rhos=[1e-3,5e-3,1e-2]
labls=[r'$\rho$'+"'="+str(x) for x in rhos]
symbs=['+-','^-','o-']
plt.subplot(1,2,1) 
for j in range(1,4): 
    KPITL(rhos[j-1],pITL)
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.xlabel('Time [s]')
plt.text(45,.6,'KP-ITL equation')
plt.ylabel('ITL [a.u]')
plt.title('(a)')
plt.subplot(1,2,2) 
Tiso=[170,180,190] 
pITL=[s*np.exp(-E/(kB*(x+273.15))) for x in Tiso]
labls=['T='+str(x)+r'$^o$C' for x in Tiso]
for j in range(1,4): 
    KPITL(0.01,pITL[j-1])
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.xlabel('Time [s]')
plt.ylabel('ITL [a.u]')
plt.title('(b)') 
plt.tight_layout()
plt.show()