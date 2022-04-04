# Deconvolution of Glocanin #1 with original MOK equation
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
from scipy import optimize
from prettytable import PrettyTable 
data = np.loadtxt('aluminaTLshort.txt')
x_data= 273+np.array(data[:, 0])
y_data=np.array(data[:, 1])/max(data[:, 1]) 
kB,  beta=  8.617e-5, 1
def TLMOK(T,A,sprime,E,alpha):
    expint=kB*(T**2.0)/(beta*E)*\
	    np.exp(-E/(kB*T))*(1-2*kB*T/E)
    Ft=np.exp(sprime*expint)
    return A*np.exp(-E/(kB*T))*Ft/((Ft-alpha)**2.0)
params, cov=optimize.curve_fit(TLMOK,x_data,y_data,([2e10,2e9,\
1.,.01]),bounds=((1e9,1e8,.9,1e-4),(1e17,1e14,1.3,1)))    
Tmax=x_data[np.argmax(y_data)]
A=format(params[0],"10.1E")
dA = format(np.sqrt(cov[0][0]),"10.1E")
sprime=format(params[1],"10.1E")
dsprime = format(np.sqrt(cov[1][1]),"10.1E")
E=round(params[2],3)
dE = round(np.sqrt(cov[2][2]),3)
alpha=round(params[3],5)
dalpha = round(np.sqrt(cov[3][3]),3)
res=TLMOK(x_data, *params)-y_data
FOM=round(100*np.sum(abs(res))/np.sum(y_data),2)
myTable = PrettyTable()
myTable.add_row([ "A","dA", "s' (s^-1)","ds' (s^-1)"]);
myTable.add_row([A,dA,sprime,dsprime]);
myTable.add_row([" "]*4);
myTable.add_row(["E (eV)","dE(eV)","alpha","dalpha"]);
myTable.add_row([E,dE,alpha,dalpha]);
print(myTable)