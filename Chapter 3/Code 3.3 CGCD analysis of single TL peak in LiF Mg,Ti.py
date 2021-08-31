# Deconvolution of LiF  TL peak with KV-TL equation
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
from prettytable import PrettyTable 
import warnings
warnings.filterwarnings("ignore")
from scipy.special import lambertw
import os
os.chdir('C:/Users/Bill Pagonis/Desktop/pythonVP')
data = np.loadtxt('LiFTL150degC.txt')
x_data,y_data = data[:, 0][0:45], data[:, 1][0:45]/max( data[:,\
1][0:45])
x_data=[273+u for u in np.array(x_data)]
x_data=np.array(x_data)
kB=8.617E-5
Imax=max(y_data)
def W_func(T,Tmax,R, E):
    F=kB*(T**2.0)*np.exp(-E/(kB*T))*(1-2*kB*T/E)/E
    Fm=kB*(Tmax**2.0)*np.exp(-E/(kB*Tmax))*(1-2*kB*Tmax/E)/E
    a=kB*Tmax**2.0*(1-1.05*R**1.26)
    Z=R/(1-R)-np.log((1-R)/R)+(F*E*np.exp(E/(kB*Tmax)))/a
    Zm=R/(1-R)-np.log((1-R)/R)+(Fm*E*np.exp(E/(kB*Tmax)))/a
    argW=np.real(lambertw(np.exp(Z)))
    argWm=np.real(lambertw(np.exp(Zm)))
    return Imax*np.exp(-E/(kB*T)*(Tmax-T)/Tmax)*\
        (argWm+argWm**2.0)/(argW+argW**2.0)  
params,cov=optimize.curve_fit(W_func,x_data,y_data,\
p0=(490,1e-6,2.0))
plt.subplot(2,1, 1)
plt.plot(x_data, W_func(x_data, *params),'-',linewidth=4)
plt.scatter(x_data, y_data, label='Experiment')
plt.plot(x_data, W_func(x_data, *params),
c='r',linewidth=3, label='KV-TL equation')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('TL signal [a.u.]')
plt.xlabel(r'Temperature T [K]')
plt.text(460, .5,'LiF:Mg,Ti')
plt.text(460, .37,'Remnant TL')
plt.subplot(2,1, 2)
plt.plot(x_data,W_func(x_data, *params)-\
y_data,'o',c='r',linewidth=2,label='Residuals')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('Residuals')
plt.xlabel(r'Temperature T [K]')
plt.ylim(-.2,.2);
plt.tight_layout()
Tmax=round(params[0],1)
R=format(params[1],"10.1E")
E=round(params[2],3)
dE = round(np.sqrt(cov[2][2]),3)
res=W_func(x_data, *params)-y_data
FOM=round(100*np.sum(abs(res))/np.sum(y_data),3)
myTable = PrettyTable([ "Tmax","R", "E(eV)","dE(eV)","FOM(%)"]) 
myTable.add_row([Tmax,R,E,dE,FOM]);
print(myTable)
plt.show()