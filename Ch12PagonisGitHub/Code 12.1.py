#Simulate TR-PL experiments with Nikiforov/Pagonis model
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
# Initial conditions.
n0, m10, m20, nc0, nF0, n3P0, n1P0 =  0, 1e14,0,0,1e14,0,0
t = np.linspace(0, 0.2, 100)
# The Pagonis model differential equations.
[E,s,alpha, delta1, delta2,Gamma,N,M1,M2,C,W,w1,w2,w3,f,kb,hr]=\
[1.3,1e13,1e-14,1e-12,1e-14,1e-11,1e13,1e15,1e14,1e13,1,1,3e3,\
29,1e10,8.617e-5,1]
def PagonisAlumina(y, t):
    n, m1, m2, nc, nF, n3P, n1P=y
    dndt= -s*np.exp(-E/(kb*(273+hr*t)))*n+alpha*(N-n)*nc
    dm1dt= delta1*(M1-m1)*nc
    dm2dt= delta2*(M2-m2)*nc
    dncdt= s*n*np.exp(-E/(kb*(273+T)))-delta1*(M1-m1)*nc-\
    delta2*(M2-m2)*nc-Gamma*nF*nc-alpha*(N-n)*nc+w1*n1P
    dnFdt= -Gamma*nF*nc+w1*n1P
    dn3Pdt=w2*n1P-C*np.exp(-W/(kb*(273+T)))*n3P-w3*n3P
    dn1Pdt=f+Gamma*nF*nc-w1*n1P-w2*n1P
    return dndt,dm1dt,dm2dt,dncdt,dnFdt,dn3Pdt,dn1Pdt
# Initial conditions vector
y0 = n0, m10, m20, nc0, nF0, n3P0, n1P0
plt.subplot(1,2, 1)
temps=[120,180,210]
cols=['r^','bo','gs']
for i in range(0,3,1):
    T=temps[i]
    ret = odeint(PagonisAlumina, y0, t)
    n, m1, m2, nc, nF, n3P, n1P = ret.T
    plt.plot(t, w3*n3P, cols[i])
plt.text(.1,8e9,'T=120$^{o}C$')
plt.text(.1,3.5e9,'T=180$^{o}C$')
plt.text(.1,1e9,'T=210$^{o}C$')
plt.text(.0,9e9,'TR-PL pulse')
plt.text(.1,4.2e9,'Stimulation')
plt.title('(a)')
plt.ylabel('TR-PL signal [a.u.]')
plt.xlabel('Time [s]')
plt.subplot(1,2, 2)
temps=range(40,240,10)
areas=[]
for i in range(len(temps)):
    T=temps[i]
    ret = odeint(PagonisAlumina, y0, t)
    n, m1, m2, nc, nF, n3P, n1P = ret.T
    u=np.sum(w3*n3P)
    areas.append(u)
plt.plot(temps,areas,'bo-')
plt.title('(b)')
plt.ylabel('Integrated TR-PL signal [a.u.]')
plt.xlabel('Stimulation temperature [$^{o}C$]')
plt.text(50,4e11,'Thermal')
plt.text(50,3.5e11,'quenching')
plt.text(50,3e11,'of intensity')
plt.tight_layout()
plt.show()