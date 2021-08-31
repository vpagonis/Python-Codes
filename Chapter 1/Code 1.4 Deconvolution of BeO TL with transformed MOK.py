# Deconvolution of BeO glow curve with transformed  MOK 
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
from prettytable import PrettyTable 
kB=8.617E-5
def MOK(T,alpha,E,Tm):
     fmok=(2.6-0.9203*alpha+0.324*(alpha**3.338))/(2.6-\
     2.9203*alpha+0.324*(alpha**3.338))
     FT=np.exp((1/fmok)*(T**2.00)/(Tm**2.0)*np.exp(E*(T-Tm)/\
     (kB*T*Tm))*(1-2.0*kB*T/E))
     FTm=np.exp((1-2*kB*Tm/E)/fmok)
     return imax*np.exp(E*(T-Tm)/(kB*T*Tm))*((FTm-alpha)**2.0)\
     	*FT/(FTm*((FT-alpha)**2.0))
import os
os.chdir('C:/Users/Bill Pagonis/Desktop/pythonVP')
data = np.loadtxt('BeOTL140.txt')
x_data= np.array(data[:, 0])+273
y_data=np.array(data[:, 1])
plt.plot(x_data,y_data,'o')  
imax=max(y_data)     #find imax from given data
T=np.arange(390,520,1);
params,cov =optimize.curve_fit(MOK,x_data,y_data,p0=(.5,1,460),
      bounds=((1e-10,.9,450),(0.99999,1.3,480)))
plt.subplot(2,1, 1)
plt.scatter(x_data, y_data, label='Data')
plt.plot(x_data, MOK(x_data,*params[0:3]),
c='r',linewidth=3, label='MOK')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('TL signal [a.u.]')
plt.xlabel(r'Temperature T [K]')
plt.text(400, 1.5e5,'BeO')
plt.subplot(2,1, 2)
plt.plot(x_data,MOK(x_data, *params[0:3])-y_data,'o',c='r',\
linewidth=2,label='Residuals')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('Residuals')
plt.xlabel(r'Temperature T [K]')
plt.tight_layout()
alphafit = format(params[0],"10.1E")
Efit = round(params[1],3)
dalpha, dE ,dTm= np.round(np.sqrt(np.diag(cov)),3)
res=MOK(x_data, *params)-y_data
FOM=round(100*np.sum(abs(res))/np.sum(y_data),2)
myTable = PrettyTable(["alpha", "E (eV)","dE (eV)",\
"FOM (%)"]) 
myTable.add_row([alphafit, Efit, dE,FOM]);
print(myTable)
plt.show()