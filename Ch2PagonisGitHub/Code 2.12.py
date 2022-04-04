# Plot of the analytical equation for GOK
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
kB,       E,  s,    beta , n0,    N=\
8.617e-5, 1, 1e12,   1 ,   1e10,  1e10
def TLgen(b,tmps):
    expint=kB*((273+tmps)**2.0)/E*\
	    np.exp(-E/(kB*(273+tmps)))*(1-2*kB*(273+tmps)/E)
    a=(n0**(b-1))*(s/(N**(b-1)))
    return n0*a*np.exp(-E/(kB*(273+tmps)))*\
    ((1+(b-1)*a/beta*expint)**(-b/(b-1)) )
tims= range(25,int(220/beta)-1,3)
temps=[beta*tim for tim in tims]
plt.plot(temps,[TLgen(1.001,x) for x in temps],'+-',
c='r',label=r'b=1.001')
plt.plot(temps,[TLgen(1.5,x) for x in temps] ,'o-',
 c='g',label=r'b=1.5')
plt.plot(temps,[TLgen(2.0,x) for x in temps] ,'^-',
 c='b',label=r'b=2')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.title('Compare GOK TL glow curves: b=1, 1.5 and 2')
plt.ylabel('TL [a.u.]')
plt.xlabel(r'Temperature T [$^{o}$C]')
plt.show()