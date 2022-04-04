# Deconvolution of Glocanin TL with the transformed FOK-TL eqt
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
from scipy import optimize
from prettytable import PrettyTable 
data = np.loadtxt('glocanin1.txt')
x_data,y_data = data[:, 0], data[:, 1]/max(data[:, 1])
kB,  beta= 8.617e-5,  1
imax=max(y_data)
Tmax=x_data[np.argmax(y_data)]
def TLFOK(T,E):      
    return imax*np.exp(1+(E/(kB*T))*((T-Tmax)/Tmax)-\
    (T**2/Tmax**2)*(1-2*kB*T/E)*np.exp((E/(kB*T))*\
    ((T-Tmax)/Tmax))-2*kB*Tmax/E)
params, cov=optimize.curve_fit(TLFOK,x_data,y_data,(1))
E=round(params[0],3)
dE = round(np.sqrt(cov[0][0]),4)
res=TLFOK(x_data, *params)-y_data
FOM=round(100*np.sum(abs(res))/np.sum(y_data),2)
myTable = PrettyTable([ "E(eV)","dE","FOM(%)"]) 
myTable.add_row([E,dE,FOM]);
print(myTable)