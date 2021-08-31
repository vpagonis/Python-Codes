# Fit TL with KP-TL (Kitis-Pagonis) analytical equation 
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
from prettytable import PrettyTable 
import os
os.chdir('C:/Users/Bill Pagonis/Desktop/pythonVP')
data = np.loadtxt('MBO6gynew.txt')
x_data,y_data = data[:, 0]+273.15, data[:, 1]
z,   kB,      En=\
1.8, 8.617E-5,1.0
Tmax=x_data[np.argmax(y_data)]
imax=max(y_data)
def test_func(T, rho):
    fm=4.90537*rho**1.21038  
    Fm=np.log(1+(1/fm)*(1-2*kB*T/En))
    F=np.log(1+(1/fm)*((T/Tmax)**2.0)*np.exp(-En*(Tmax-T)/(kB*T*Tmax))*(1-2*kB*T/En))
    return imax*np.exp(-En*(Tmax-T)/(kB*T*Tmax))*((F/Fm)**2.0)*\
    np.exp(-rho*(F**3.0)-F)/np.exp(-rho*(Fm**3.0)-Fm)
params, cov = optimize.curve_fit(test_func,\
x_data, y_data,bounds=(0.002,.02))
drho= round(np.sqrt(cov[0][0]),5)
plt.subplot(2,1, 1)
plt.plot(x_data, y_data,'o', c='lightgreen',label='Experiment')
plt.plot(x_data, test_func(x_data, *params),label='KP-TL equation',
         c='black',linewidth=2)
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
plt.ylim(-1,1);
plt.ylabel('Residuals')
plt.xlabel(r'Temperature T [$^{o}$C]')
plt.tight_layout()
rho=round(params[0],4)
res=test_func(x_data, *params)-y_data
FOM=round(100*np.sum(abs(res))/np.sum(y_data),2)
myTable=PrettyTable([ "rho",  "d(rho)","FOM"]);  
myTable.add_row([rho,drho, FOM]);
#print(myTable)
plt.show()