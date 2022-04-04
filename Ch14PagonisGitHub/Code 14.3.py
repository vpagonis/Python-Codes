#plot dose response  PKC-S equation for different beta, B
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import lambertw
import warnings
warnings.filterwarnings("ignore")   
def PKCS(B,beta):       
    u=np.real(lambertw(B*np.exp(np.abs(B)-(D/Dc))))/B
    u=(1-(u**beta))
    plt.plot(D,u,symbs[j-1], linewidth=2,
    label=labls[j-1])
    plt.yscale("log")
    plt.xscale("log")
B, Dc= 2, 2
D = np.linspace(.1, 40, 100)
betas=[.01,.1,.9]
labls=[r'$\beta$='+str(x) for x in betas]
symbs=['+-','^-','o-']
plt.subplot(1,2,1) 
for j in range(1,4): 
    PKCS(B,betas[j-1])
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.xlabel('Dose [Gy]')
plt.text(.1,.8,'PKC-S equation')
# plt.text(45,.6,'PKC equation')
plt.ylabel('n(D)/N')
plt.title('(a)')
plt.subplot(1,2,2) 
Bs=[.1,1,2]
labls=['B='+str(x) for x in Bs]
for j in range(1,4): 
    PKCS(Bs[j-1],.01)
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.xlabel('Dose [Gy]')
plt.ylabel('n(D)/N')
plt.title('(b)') 
plt.tight_layout()
plt.show()