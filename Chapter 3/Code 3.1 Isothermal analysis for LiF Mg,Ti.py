# Isothermal analysis for LiF:Mg,Ti
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from prettytable import PrettyTable 
from scipy import stats
def expon(x,A,tau):    
    return A*np.exp(-x/tau)
import os
os.chdir('C:/Users/Bill Pagonis/Desktop/pythonVP')
files=('LiF180.txt','LiF190.txt','LiF210.txt','LiF220.txt')
slopes=[0 for x in range(4)]
plt.subplot(1,2, 1)
plt.title('(a')
plt.ylim(0,1.2);
plt.ylabel('Isothermal signal [a.u.]')
plt.xlabel(r'Time t [s]')
sym=(['o','s','^','x'])
lab=('180$^{o}$C','190$^{o}$C','200$^{o}$C','210$^{o}$C')
for i in range(4):
    dat = np.loadtxt(files[i])
    x_data,y_data = dat[:, 0],  dat[:, 1]
    plt.plot(x_data, y_data,sym[i],label=lab[i])
    params, cov=optimize.curve_fit(expon,x_data,y_data)
    print(params)
    slopes[i]=1/params[1]
    plt.plot(x_data, expon(x_data, *params))
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
temps=([180,190,200,210])
xArrh=[1/(8.617e-5*(x+273)) for x in temps]
yArrh=np.log(slopes)
slope, intercept,r_value,p_value, std_err=\
stats.linregress(xArrh,yArrh)
plt.subplot(1,2, 2)
plt.plot(xArrh,yArrh,'o')
fit=[x*slope+intercept for x in xArrh]
plt.plot(xArrh,fit,label='Fit')
plt.xlabel('1/(kT) (eV$^{-1}$)')
plt.ylabel(r'ln[$\lambda$(T)]')
plt.text(24.5,-3,'E=-slope')
plt.title('(b)');
plt.tight_layout()
slope, intercept, rsqr, p_value,std_err=\
np.round(slope,2),np.round(intercept,2),np.round(\
r_value**2.0,3),format(p_value,"10.2E"), np.round(std_err,2)
myTable=PrettyTable(['E (eV)','dE (eV)','p-value','R^2']) 
myTable.add_row([-slope, std_err, p_value,rsqr]);
print(myTable)
plt.show()