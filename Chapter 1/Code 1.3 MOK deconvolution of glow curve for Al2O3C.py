# Deconvolution of Al2O3:C peak with MOK 
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from prettytable import PrettyTable 
import os
os.chdir('C:/Users/Bill Pagonis/Desktop/pythonVP')
data = np.loadtxt('aluminaTL2.txt')
x_data= np.array(data[:, 0])
y_data=np.array(data[:, 1])/max(data[:, 1]) 
x_data=[273+u for u in np.array(x_data)]
x_data=np.array(x_data)
plt.plot(x_data,y_data,'o')
kB,  beta=  8.617e-5, 1
def TLMOK(T,A,sprime,E,alpha):
    expint=kB*(T**2.0)/E*\
	    np.exp(-E/(kB*T))*(1-2*kB*T/E)
    Ft=np.exp(sprime*expint)
    return A*np.exp(-E/(kB*T))*Ft/((Ft-alpha)**2.0)
params, cov=optimize.curve_fit(TLMOK,x_data,y_data,([2e10,2e9,\
1.,.01]),bounds=((1e9,1e8,.9,1e-4),(1e17,1e14,1.3,1)))    
plt.subplot(2,1,1)
plt.plot(x_data, TLMOK(x_data, *params),'-',linewidth=4)
plt.scatter(x_data, y_data, label='Experiment')
plt.plot(x_data, TLMOK(x_data, *params),
c='r',linewidth=3, label='MOK equation')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('TL signal [a.u.]')
plt.xlabel(r'Temperature T [K]')
plt.text(520, .7,r'Al$_{2}$O$_{3}$:C')
plt.subplot(2,1,2)
plt.plot(x_data,TLMOK(x_data, *params)-\
y_data,'o',linewidth=2,label='Residuals')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('Residuals')
plt.xlabel(r'Temperature T [K]')
plt.ylim(-.1,.1);
plt.plot(([min(x_data),max(x_data)]),([0,0]),'r-')
plt.tight_layout()
Tmax=x_data[np.argmax(y_data)]
A=format(params[0],"10.1E")
dA = format(np.sqrt(cov[0][0]),"10.1E")
sprime=format(params[1],"10.1E")
E=round(params[2],3)
dE = round(np.sqrt(cov[2][2]),3)
alpha=round(params[3],5)
res=TLMOK(x_data, *params)-y_data
FOM=round(100*np.sum(abs(res))/np.sum(y_data),2)
myTable = PrettyTable([ "A", "s' (K^-1)","E(eV)","dE(eV)",\
"alpha","FOM(%)"]) 
myTable.add_row([A,sprime,E,dE,alpha,FOM]);
print(myTable)
plt.show()