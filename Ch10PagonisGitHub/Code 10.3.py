#Example of summing partial LM-IRSL curves in the EST model
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
##### function to partial CW curves ###
def partialLMIRSL(tims,rprime):
    seff=A*np.exp(-(rho**(-1/3))*rprime)
    LMIRSLpartial=3*(rprime**2.0)*np.exp(-(rprime**3.0))*seff*\
    tims*np.exp(-seff*tims**2/(2*P))
    return LMIRSLpartial*dr
A, P=5, 200
t = np.linspace(0, P, P)
rho=0.013               		# rho-prime value
dr=0.1
rprimes=np.arange(0,2.2,dr)    # rprime=0-2.2
plt.subplot(1,2, 1)
for i in range(len(rprimes)):
    plt.plot(t,[partialLMIRSL(x,rprimes[i]) for x in t])
plt.text(80,.35,"Partial LM-IRSL")
plt.text(80,.30,"for r'=0-2.2")
plt.ylabel('Partial LM-IRSL signal [a.u.]')
plt.xlabel(r'Time [s]')
plt.title('(a)')
plt.text(80,.25,'EST model')
plt.subplot(1,2,2)
u=np.array([[partialLMIRSL(x,rprimes[i]) for x in t]\
for i in range(len(rprimes))])
plt.plot(t,sum(u),c='b')
plt.text(80,1.8, 'Sum of')
plt.text(80,1.55, 'partial')
plt.text(80,1.3, 'LM-IRSL curves')
plt.ylabel('LM-IRSL signal [a.u.]')
plt.xlabel(r'Time [s]')
plt.title('(b)')
plt.tight_layout()
plt.show()