#Deconvolution of 5-peak glow curve for BAL21 sample
# using the  KP-TL equation
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
from prettytable import PrettyTable 
import warnings
warnings.filterwarnings("ignore")
import os
os.chdir('C:/Users/Bill Pagonis/Desktop/pythonVP')
data = np.loadtxt('initialsensBAL21.txt')
x_data,y_data = data[:, 0], data[:, 1]
plt.scatter(x_data, y_data, label='Experiment')
kB=8.617E-5
nPks=5
z,   kB,     s, =1.8, 8.617E-5,1e12
def TL(T, B,En ,rho):
    return abs(B)* np.exp(-rho*( (np.log(1+z*s*kB*(((T+\
			273)**2.0)/np.abs(En))*np.exp(-En/(kB*(T+273)))*\
			(1-2*kB*(T+273)/En)))**3.0))*(En**2.0-6*(kB**2.0)*\
    	((T+273)**2.0))*((np.log(1+z*s*kB*(((T+273)**2.0)/\
			abs(En))*np.exp(-En/(kB*(T+273)))*(1-2*kB*(T+273)/\
			En)))**2.0)/(En*kB*s*((T+273)**2)*z-2*(kB**2.0)*\
			s*z*((T+273)**3.0)+np.exp(En/(kB*(T+273)))*En)
def total_TL(T, *inis): 
    u=np.array([0 for i in range(len(x_data))])
    Bs, rho=    inis[0:nPks],inis[-1]
    for i in range(nPks):
        u=u+TL(T,Bs[i],Ens[i],rho)
    return u
inis=(5e16,9e16,10e16,1e16,.1e16,.004)
Ens=[.82, .95,1.06,1.19,1.355]

params, params_covariance = optimize.curve_fit(total_TL,\
x_data,y_data,p0=inis)
plt.subplot(2,1, 1)
plt.scatter(x_data, y_data, label='BAL21 feldspar data')
plt.plot(x_data, total_TL(x_data, 
*params),c='r',linewidth=3, label='Analytical KV-TL')
plt.plot(x_data, TL(x_data, params[0],Ens[0],params[5]))
for i in range(0,nPks):
     plt.plot(x_data, TL(x_data, params[i],Ens[i],params[5]))
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('TL signal [a.u.]')
plt.xlabel(r'Temperature T [K]')
plt.subplot(2,1, 2)
plt.scatter(x_data,total_TL(x_data, *params)
   -y_data,c='r',linewidth=2,label='Residuals')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('Residuals')
plt.xlabel(r'Temperature T [K]')
plt.ylim(-5e4,5e4);
plt.tight_layout()
res=total_TL(x_data, *params)-y_data
FOM=round(100*np.sum(abs(res))/np.sum(y_data),2)

print('FOM=',FOM,' %')
plt.show()