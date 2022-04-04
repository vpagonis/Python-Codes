#plot coefficient f(D) for different beta, B
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import lambertw
import warnings
warnings.filterwarnings("ignore")   
def PKCS(B,beta):       
    u=np.real(lambertw(B*np.exp(np.abs(B)-(D/Dc))))/B
    u=(1-(u**beta))/D
    u1=np.real(lambertw(B*np.exp(np.abs(B)-(.1/Dc))))/B    
    u1=(1-(u1**beta))/.1
    plt.plot(D,u/u1,symbs[j-1], linewidth=2,
    label=labls[j-1])
    plt.xscale("log")
B, Dc= 2, 2
D = np.linspace(.1, 40, 40)
betas=[.1,.3,.9]
labls=[r'$\beta$='+str(x) for x in betas]
symbs=['+-','^-','o-']
plt.subplot(1,2,1) 
for j in range(1,4): 
    PKCS(B,betas[j-1])
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.xlabel('Dose [Gy]')
plt.ylabel('f(D)')
plt.title('(a)')
plt.subplot(1,2,2) 
Bs=[1,2,3]
labls=['B='+str(x) for x in Bs]
for j in range(1,4): 
    PKCS(Bs[j-1],.1)
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylim(0,3)
plt.xlabel('Dose [Gy]')
plt.ylabel('f(D)')
plt.title('(b)') 
plt.tight_layout()
plt.show()