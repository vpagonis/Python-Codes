# Deconvolution of TL user data (.txt file, GOK)
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
from prettytable import PrettyTable
import warnings
warnings.filterwarnings("ignore")
data = np.loadtxt('lbodata.txt')
x_data,y_data = data[:, 0], data[:, 1]
kB=8.617E-5
def TL(T, imax,b, En,Tmax):
    return imax* np.exp(En/(kB*T)*(T-Tmax)/Tmax)*(b**\
		    ((b/(b-1))))*((1+(b-1)*2*kB*Tmax/En+(b-1)*(1-2*kB*T/En)*\
		    np.exp(En/(kB*T)*(T-Tmax)/Tmax)*(T**2.0)/(Tmax**2.0))**\
		    (b/(1-b)))
def total_TL(T, imax1,b1, En1,Tmax1, imax2,b2, En2,Tmax2):
    return TL(T, imax1,b1, En1,Tmax1)+TL(T,imax2,b2, En2,Tmax2)                           
params, cov = optimize.curve_fit(total_TL, x_data,
y_data,bounds=((80,1.001,.8,430, 20,1.001,.8,480),
(140,2.0,1.3,480, 40,2.0,1.3,540)))
imax1,b1,E1,Tmax1=int(params[0]),round(params[1],2),\
round(params[2],2),int(params[3]),
imax2,b2,E2,Tmax2=int(params[4]),round(params[5],2),\
round(params[6],2),int(params[7]),
beta= 1
s1=np.exp(E1/(kB*Tmax1))*(beta*E1/(kB*(Tmax1**2.0)))\
/(1+(b1-1)*2*kB*Tmax1/E1)
s1=format(s1,"10.2E")
s2=np.exp(E2/(kB*Tmax2))*(beta*E2/(kB*(Tmax2**2.0)))\
/(1+(b2-1)*2*kB*Tmax2/E2)
s2=format(s2,"10.2E")
db1 = round(np.sqrt(cov[1][1]),2)
dE1 = round(np.sqrt(cov[2][2]),2)
db2 = round(np.sqrt(cov[5][5]),2)
dE2 = round(np.sqrt(cov[6][6]),2)
res=total_TL(x_data, *params)-y_data
FOM=round(100*np.sum(abs(res))/np.sum(y_data),2)
myTable = PrettyTable(["Imax", "b","db", "E(eV)",\
"dE(eV)", "s(s^-1)","FOM(%)"]) 
myTable.add_row([imax1,b1,db1,E1,dE1,s1,FOM]);
myTable.add_row([imax2,b2,db2,E2,dE2,s2,"-"]);
print(myTable)