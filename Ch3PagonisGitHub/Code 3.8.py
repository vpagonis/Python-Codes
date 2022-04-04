# Deconvolution of Glocanin TL with transformed  MOK 
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
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
data = np.loadtxt('aluminaTLshort.txt')
x_data= 273+np.array(data[:, 0])
y_data=np.array(data[:, 1])
y_data=y_data/max(y_data)
imax=max(y_data)     #find imax from given data
Tmax=x_data[np.argmax(y_data)]
T=np.arange(390,570,1)
params,cov =optimize.curve_fit(MOK,x_data,y_data,p0=(.9,1.07,\
Tmax),bounds=((1e-5,.9,Tmax-5),(0.99999,1.3,Tmax+5)))
alphafit,Efit,Taxfit = np.round(params,2)
dalpha, dE,dTm = np.round(np.sqrt(np.diag(cov)),2)
res=MOK(x_data, *params)-y_data
FOM=round(100*np.sum(abs(res))/np.sum(y_data),2)
myTable=PrettyTable(["alpha","daplha","E (eV)","dE","FOM(%)"]) 
myTable.add_row([alphafit,dalpha, Efit, dE,FOM])
print(myTable)