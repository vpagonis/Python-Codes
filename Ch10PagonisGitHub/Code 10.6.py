#Fitting CW-IRSL signal with stretched exponential
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
# pLM in this code represents the effective frequency s'  
def KPLM(tau,beta):    
    SE=A*np.exp(-(t/tau)**beta)   
    plt.plot(t,SE/max(SE),symbs[j-1], linewidth=2,
    label=labls[j-1])
t = np.linspace(1, 100, 30)
A=1
taus=[5,10,15]
beta=0.5
labls=[r'$\tau$'+"="+str(x)+" s" for x in taus]
symbs=['+-','^-','o-']
plt.subplot(1,2,1) 
for j in range(1,4): 
    KPLM(taus[j-1],beta)
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.xlabel('Time [s]')
plt.text(50,.6,'Stretched')
plt.text(50,.52,'Exponential')
plt.ylabel('LM-IRSL [a.u]')
plt.title('(a)')
plt.subplot(1,2,2) 
betas=[.3,.5,.7]
tau=15
labls=[r"$\beta$="+str(x) for x in betas]
for j in range(1,4): 
    KPLM(tau,betas[j-1])
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.xlabel('Time [s]')
plt.ylabel('LM-IRSL [a.u]')
plt.title('(b)') 
plt.tight_layout()
plt.show()