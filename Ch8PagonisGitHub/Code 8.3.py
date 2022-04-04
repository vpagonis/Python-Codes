# Plots  of the Bulur FOK equation for LM-OSL
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
def BulurFOK(P,lamda,j):
    plt.plot(t,n0*lamda*t/P*np.exp(-lamda*t**2.0/(2*P)),\
    symbs[j],label=labls[j])
n0 = 1e10
t = np.linspace(0, 100, 100)
plt.subplot(1,2, 1)
symbs=['o-','+-','^-']
labls=[r'$\lambda$=0.1 s$^{-1}$',r'$\lambda$=0.2 s$^{-1}$',\
r'$\lambda$=0.3 s$^{-1}$']
BulurFOK(100,0.1,0)
BulurFOK(100,0.2,1)
BulurFOK(100,0.3,2)
plt.title('(a)')
plt.text(60,1.9e8,'Bulur eqt.')
plt.text(60,1.7e8,'for LM-OSL')
plt.text(60,1.5e8,r'Variable $\lambda$')
plt.xlabel('Time [s]')
plt.ylabel('CW-OSL [a.u.]')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.subplot(1,2, 2)
labls=['P=200 s','P=300 s','P=400 s'] 
BulurFOK(200,.3,0)
BulurFOK(300,.3,1)
BulurFOK(400,.3,2)
plt.title('(b)')
plt.text(60,1.6e8,'Variable P')
plt.xlabel('Time [s]')
plt.ylabel('CW-OSL [a.u.]')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.tight_layout()
plt.show()