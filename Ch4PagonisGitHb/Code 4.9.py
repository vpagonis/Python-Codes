# Dose response n(t) in the TA-EST model by Brown et al.
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
##### function to find distribution of distances ###
def findDistr(tim):
    P=Po*np.exp(-(rho**(-1/3))*rprimes)
    Peff=(P*s)*np.exp(-E/(kB*(273+Tirr)))/(P+s)
    return 3*(rprimes**2.0)*np.exp(-(rprimes**3.0))*\
    X/(Do*Peff+X)*(1-np.exp(-(X/Do+Peff)*tim))
Tirr=-4                   # irradiation temeprature
Po=2e15                   # tunneling frequency factor
s=2e15                    # thermal frequency factor
rho=1e-3               	  # rho-prime values 0.005-0.02
Do=1600                   # D0 in Gy
yr=365*3600*24            # year in seconds
X=2.85/(1000*yr)          # natural dose rate 2.85 Gy/ka
times=[yr*1e4,yr*1e5,yr*.2e6,yr*.4e6,yr*.7e6,yr*1e6] # time (s)
t=[times[x]/(1e6*yr) for x in range(6)]
kB, E =8.617e-5, 1.3
dr=0.05
rprimes=np.arange(0,2.2,dr)    # rprime=0-2.2
distribs=[findDistr(x) for x in times]
plt.subplot(1,2, 1)
labls=[str(x) for x in t]
for j in range(6): 
    plt.plot(rprimes,distribs[j],'-',label=labls[j])
plt.text(1.55,.32,'Irradiation')
plt.text(1.55,.27,'time in Ma')
plt.xlim(0,2.4)
plt.ylabel('Trapped electrons distribution n(r$^{\'}$,t)/N')
plt.xlabel('Distance r$^{\'}$')
plt.title('(a)')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.subplot(1,2, 2)
sums=[sum(distribs[x])*dr for x in range(6)]
plt.plot(t,sums,'o-')
plt.title('(b)')
plt.ylabel('Total trapped electrons n(t)/N')
plt.xlabel('Irradiation time [Ma]')
plt.text(.33,.16,'TA-EST model')
plt.tight_layout()
plt.show()