# Deconvolution of Al2O3:C glow curve (GOK)
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
from prettytable import PrettyTable 
import warnings
warnings.filterwarnings("ignore")
import os
os.chdir('C:/Users/Bill Pagonis/Desktop/pythonVP')
data = np.loadtxt('aluminaTL2.txt')
x_data,y_data = data[:, 0], data[:, 1]
x_data=[273+u for u in np.array(x_data)]
x_data=np.array(x_data)
kB=8.617E-5
imax=max(y_data)
def GOK_func(T, Tmax,b, En):
    return imax* np.exp(En/(kB*T)*(T-Tmax)/Tmax)*(b**\
    ((b/(b-1))))*((1+(b-1)*2*kB*Tmax/En+(b-1)*(1-2*kB*T/\
    En)*np.exp(En/(kB*T)*(T-Tmax)/Tmax)*(T**2.0)/(Tmax**\
    2.0))**(b/(1-b)))                                   
params, cov=optimize.curve_fit(GOK_func,x_data,\
y_data,bounds=((460,1.001,1),(520,2.0,1.3)))
plt.subplot(2,1, 1)
plt.scatter(x_data, y_data, label='Experiment')
plt.plot(x_data, GOK_func(x_data, *params),
c='r',linewidth=3, label='GOK equation')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('TL signal [a.u.]')
plt.xlabel(r'Temperature T [K]')
plt.text(400, 1e6,'Al$_{2}$O$_{3}$:C')
plt.xlim(375,550);
plt.subplot(2,1, 2)
plt.plot(x_data,GOK_func(x_data, *params)-\
y_data,'o',c='r',linewidth=2,label='Residuals')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('Residuals')
plt.xlabel(r'Temperature T [K]')
plt.ylim(-4e5,4e5,'o');
plt.xlim(375,550);
plt.tight_layout()
Tmax=params[0]
s=np.exp(params[2]/(kB*Tmax))*(params[2]/(kB*(Tmax**2.0)))
sf=format(s, "10.1E")
Tmax=format(params[0],"10.1E")
b=round(params[1],3)
E=round(params[2],4)
db = round(np.sqrt(cov[1][1]),3)
dE = round(np.sqrt(cov[2][2]),3)
res=GOK_func(x_data, *params)-y_data
FOM=round(100*np.sum(abs(res))/np.sum(y_data),2)
myTable = PrettyTable([ "b","db", "E(eV)","dE(eV)",\
"s(s^-1)","FOM(%)"]) 
myTable.add_row([b,db,E,dE,sf,FOM]);
print(myTable)
plt.show()