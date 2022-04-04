#Example of summing partial CW-IRSL curves in the EST model
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
##### function to find distribution of distances ###
def partialCWIRSL(tims,rprime):
    seff=A*np.exp(-(rho**(-1/3))*rprime)
    CWIRSLpartial=3*(rprime**2.0)*np.exp(-(rprime**3.0))*seff*\
    np.exp(-seff*tims)
    return CWIRSLpartial*dr
A, P=5, 400
t = np.linspace(0, P, P)
rho=1e-3               		# rho-prime value
dr=0.1
rprimes=np.arange(0,2.2,dr)    # rprime=0-2.2
plt.subplot(1,2, 1)
for i in range(len(rprimes)):
    plt.plot(t,[partialCWIRSL(x,rprimes[i]) for x in t])
plt.text(120,.0035,"Partial CW-IRSL")
plt.text(120,.003,"for r'=0-2.2")
plt.ylabel('Partial CW-IRSL signal [a.u.]')
plt.xlabel(r'Time [s]')
plt.title('(a)')
plt.text(120,.005,'EST model')
plt.subplot(1,2,2)
u=np.array([[partialCWIRSL(x,rprimes[i]) for x in t]\
for i in range(len(rprimes))])
plt.plot(t,sum(u),c='b')
plt.text(100,.015, 'Sum of')
plt.text(100,.013, 'partial')
plt.text(100,.011, 'CW-IRSL curves')
plt.ylabel('CW-IRSL signal [a.u.]')
plt.xlabel(r'Time [s]')
plt.title('(b)')
plt.tight_layout()
plt.show()