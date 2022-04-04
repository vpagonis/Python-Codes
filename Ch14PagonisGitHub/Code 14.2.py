#plot dose response using PKC equation for different R, Dc
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import lambertw
import warnings
warnings.filterwarnings("ignore")   
def PKC(R,Dc):       
    u=(1+np.real(lambertw((R-1)*np.exp((R-1)-\
    (D/np.abs(Dc))))/(1-R)))  
    plt.plot(D,u/max(u),symbs[j-1], linewidth=2,
    label=labls[j-1])
Dc=20
D = np.linspace(0, 100, 100)
Rs=[.01,.5,.9]
labls=["R="+str(x) for x in Rs]
symbs=['+-','^-','o-']
plt.subplot(1,2,1) 
for j in range(1,4): 
    PKC(Rs[j-1],Dc)
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.xlabel('Dose [Gy]')
plt.text(45,.6,'PKC equation')
plt.ylabel('n(D)/N')
plt.title('(a)')
plt.subplot(1,2,2) 
Dcs=[5,10,20]
labls=['Dc='+str(x)+' Gy' for x in Dcs]
for j in range(1,4): 
    PKC(0.01,Dcs[j-1])
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.xlabel('Dose [Gy]')
plt.ylabel('n(D)/N')
plt.title('(b)') 
plt.tight_layout()
plt.show()