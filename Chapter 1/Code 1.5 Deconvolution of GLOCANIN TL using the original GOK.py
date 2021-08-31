# Deconvolution of GLOCANIN TL using the original GOK
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from prettytable import PrettyTable 
import os
os.chdir('C:/Users/Bill Pagonis/Desktop/pythonVP')
data = np.loadtxt('glocanin1.txt')
x_data,y_data = data[:, 0], data[:, 1]/max(data[:, 1])
kB= 8.617e-5
def TLGOK(T,A,sprime,E,b):    
    expint=kB*(T**2.0)/E*\
	    np.exp(-E/(kB*T))*(1-2*kB*T)/E
    return A*np.exp(-E/(kB*T))*((1+sprime*expint)**(-b/(b-1)) )
params, cov=optimize.curve_fit(TLGOK,x_data,y_data,([1,2e10,1.,\
1.5]),bounds=((1e-10,1e8,1,1),(1e20,1e14,1.3,1.999)))
plt.subplot(2,1, 1)
plt.scatter(x_data, y_data, label='Experiment')
plt.plot(x_data, TLGOK(x_data, *params),
c='r',linewidth=3, label='GOK equation')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('TL signal [a.u.]')
plt.xlabel(r'Temperature T [K]')
plt.text(400, .5,'Glocanin')
plt.text(400, .37,'Project')
plt.xlim(375,550);
plt.subplot(2,1, 2)
plt.plot(x_data,TLGOK(x_data, *params)-\
y_data,c='r',linewidth=2,label='Residuals')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('Residuals')
plt.xlabel(r'Temperature T [K]')
plt.ylim(-1e-3,1e-3);
plt.xlim(375,550);
plt.tight_layout()
A=format(params[0],"10.1E")
dA = format(np.sqrt(cov[0][0]),"10.1E")
sprime=format(params[1],"10.1E")
E=round(params[2],3)
dE = round(np.sqrt(cov[2][2]),7)
b=round(params[3],3)
db = round(np.sqrt(cov[3][3]),7)
res=TLGOK(x_data, *params)-y_data
FOM=round(100*np.sum(abs(res))/np.sum(y_data),2)
myTable = PrettyTable([ "A","s' (K^-1)","E(eV)",\
"b","FOM(%)"]) 
myTable.add_row([A,sprime,E,b,FOM]);
print(myTable)
plt.show()