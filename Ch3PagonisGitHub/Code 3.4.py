# Deconvolution of Glocanin TL #1 with transformed GOK
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
from prettytable import PrettyTable 
import warnings
warnings.filterwarnings("ignore")
data = np.loadtxt('glocanin1.txt')
x_data,y_data = data[:, 0], data[:, 1]
y_data=y_data/max(y_data)
kB=8.617E-5
imax=max(y_data)
def GOK_func(T, Tmax,b, En):
    return imax* np.exp(En/(kB*T)*(T-Tmax)/Tmax)*(b**\
    ((b/(b-1))))*((1+(b-1)*2*kB*Tmax/En+(b-1)*(1-2*kB*T/\
    En)*np.exp(En/(kB*T)*(T-Tmax)/Tmax)*(T**2.0)/(Tmax**\
    2.0))**(b/(1-b)))                                   
params, cov=optimize.curve_fit(GOK_func,x_data,\
y_data,bounds=((460,1.001,.7),(520,2.0,1.3)))
Tmax=params[0]
s=np.exp(params[2]/(kB*Tmax))*(params[2]/(kB*(Tmax**2.0)))
sf=format(s, "10.1E")
Tmax=format(params[0],"10.1E")
b=round(params[1],3)
E=round(params[2],4)
db = format(np.sqrt(cov[1][1]),"10.1E")
dE = format(np.sqrt(cov[2][2]),"10.1E")
res=GOK_func(x_data, *params)-y_data
FOM=round(100*np.sum(np.abs(res))/np.sum(y_data),2)
myTable = PrettyTable([ "b","db", "E(eV)","dE(eV)",\
"s(s^-1)"]) 
myTable.add_row([b,db,E,dE,sf]);
print("FOM=",FOM," %")
print(myTable)