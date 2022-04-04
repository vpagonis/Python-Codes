# Simultaneous irradiation and anomalous fading in nature
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
s=3e15                  # frequency factor
Do=538                  # D0 in Gy
X=3/(1000*365*3600*24)  # natural dose rate 3 Gy/ka
years=1000              # number of years elapsed
elapsedt=3.154e7*years  # change years to seconds
times=np.arange(1,elapsedt,1e6)          # times in seconds
rhos=[1e-6,2e-6,3e-6]                    #rho-prime values
##### function to find n(t)
def n(rho):
    return (1-np.exp(-dose/Do))*np.exp(-rho*(np.log(Do*s/X)\
	    **3.0))
dose=np.arange(1,3500,200)
ns=[n(x) for x in rhos] 
plt.plot(dose,ns[0],'-',c='black',
linewidth=2,label=r'$\rho$'+"'"+'=1x10$^{-6}$')
plt.plot(dose,ns[1],'-.',c='r',
linewidth=2,label=r'$\rho$'+"'"+'=2x10$^{-6}$')
plt.plot(dose,ns[2],'--',c='b',
linewidth=2,label=r'$\rho$'+"'"+'=3x10$^{-6}$')
leg = plt.legend()
plt.ylim(0,1.20)
plt.text(1500,1.1,'Acceptor density')
leg.get_frame().set_linewidth(0.0)
plt.xlabel('Natural dose [Gy]')
plt.ylabel('L$_{FADED}$(D)')
plt.title('Simultaneous irradiation and anomalous fading')
plt.show() 