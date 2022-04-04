# Fit TL with KP-TL (Kitis-Pagonis) analytical equation 
from scipy import optimize
from sympy import *
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
from prettytable import PrettyTable 
data = np.loadtxt('MBO6gynew.txt')
x_data,y_data = data[:, 0], data[:, 1]
z,   kB,      En=\
1.8, 8.617E-5,1.0
def test_func(x, B,rho, s):
    return B* np.exp(-rho*( (np.log(1+z*s*kB*(((x+\
		    273)**2.0)/np.abs(En))*np.exp(-En/(kB*(x+273)))*\
		    (1-2*kB*(x+273)/En)))**3.0))*(En**2.0-6*(kB**2.0)*\
    ((x+273)**2.0))*((np.log(1+z*s*kB*(((x+273)**2.0)/\
		    abs(En))*np.exp(-En/(kB*(x+273)))*(1-2*kB*(x+273)/\
		    En)))**2.0)/(En*kB*s*((x+273)**2)*z-2*(kB**2.0)*\
		    s*z*((x+273)**3.0)+np.exp(En/(kB*(x+273)))*En)
params, cov = optimize.curve_fit(test_func,\
x_data, y_data,bounds=(0,[1e20,.02,1e14]))
drho= round(np.sqrt(cov[1][1]),5)
plt.subplot(2,1, 1)
plt.plot(x_data, y_data,'o', c='lightgreen',label='Experiment')
plt.plot(x_data, test_func(x_data, *params),
label='KP-TL equation',   c='black',linewidth=2)
plt.title('(a)')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('TL signal [a.u.]')
plt.xlabel(r'Temperature T [$^{o}$C]')
#plt.text(400, 1.8e4,'KST4 feldspar')
plt.subplot(2,1, 2)
plt.plot(x_data,test_func(x_data, *params)-y_data,'o',\
label='Residuals')
plt.title('(b)')
plt.ylim(-1,1)
plt.ylabel('Residuals')
plt.xlabel(r'Temperature T [$^{o}$C]')
plt.tight_layout()
B,rho, s=format(params[0],"10.2E"),round(params[1],5),\
format(round(params[2],2),"10.2E")
res=test_func(x_data, *params)-y_data
FOM=round(100*np.sum(abs(res))/np.sum(y_data),2)
myTable=PrettyTable(["B", "rho",  "d(rho)",\
"s(s^-1)","FOM"])  
myTable.add_row([B,rho,drho, s,FOM])
print(myTable)
plt.show()