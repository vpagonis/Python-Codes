# Apply the heating rate method to find E,s 
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
from scipy import stats
from prettytable import PrettyTable 
kB,       E, s,     n0 =\
8.617e-5, 1, 1e12, .2e10
def TLfirst(beta,temps):
    return n0*(s/beta)*np.exp(-E/(kB*(273+temps)))*\
    np.exp(-s*kB*((273+temps)**2.0)/(beta*E)*\
    np.exp(-E/(kB*(273+temps)))*(1-2*kB*(273+temps)/E))
temps=np.arange(50,160)
It1=[TLfirst(1,x) for x in temps]
It2=[TLfirst(2,x) for x in temps]
It3=[TLfirst(3,x) for x in temps]
It4=[TLfirst(4,x) for x in temps]
plt.subplot(1,2, 1)
symb=('o-','^-','+-','x-')
for j in range(1,5,1):
    plt.plot(temps,[TLfirst(j,x) for x in temps],symb[j-1])
plt.ylabel('TL signal [a.u.]')
plt.xlabel(r'Temperature T [$^{o}$C]')
plt.text(130,6e7,r'$\beta$=1-4 K/s')
plt.xlim(70,160)
plt.ylim(0,7e7)
plt.title('(a)')
Tmax=[273+temps[np.argmax(It1)],273+temps[np.argmax(It2)],
273+temps[np.argmax(It3)],273+temps[np.argmax(It4)]]
y=[np.log(Tmax[0]**2/1),np.log(Tmax[1]**2/2),np.log(Tmax[2]**2/\
3),np.log(Tmax[3]**2/4)]
x=[1.0/(kB*u) for u in Tmax]
plt.subplot(1,2, 2)
plt.plot(x,y,"o",c="r")
plt.ylabel(r'ln(T$_{m}^{2}$/$\beta$)')
plt.xlabel(r'1/(kT$_{m}$) [eV$^{-1}$]')
x=np.array(x)
y=np.array(y)
slope, intercept,r_value,p_value,std_err=stats.linregress(x,y)
plt.plot(x,x*slope+intercept,label='Fit')
slope, intercept, rsqr, std_err=np.round(slope,2),\
np.round(intercept,2),np.round(r_value**2.0,3),\
np.round(std_err,2)
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.title('(b)')
plt.text(29,11.6,'Heating rate')
plt.text(29,11.5,'method')
plt.tight_layout()
s=format(np.exp(-intercept)*slope/kB, "10.2E")
myTable=PrettyTable([ "E (eV)",  "dE (eV)",\
"intercept I","s (s^-1)","R^2"])  
myTable.add_row([slope, std_err, intercept, s, rsqr])
print(myTable) 
plt.show()