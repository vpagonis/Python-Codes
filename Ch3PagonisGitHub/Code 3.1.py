# Deconvolution of Glocanin TL with the original FOK-TL equation
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
from scipy import optimize
from prettytable import PrettyTable 
data = np.loadtxt('glocanin1.txt')
x_data,y_data = data[:, 0], data[:, 1]/max(data[:, 1])
kB,  beta= 8.617e-5,  1
# function for evaluating the FOK-TL (R-W) equation
def TLFOK(T,A,s,E):    
    expint=kB*(T**2.0)/(beta*E)*\
	    np.exp(-E/(kB*T))*(1-2*kB*T)/E
    return A*np.exp(-E/(kB*T))*np.exp(-(s/beta)*expint) 
inis=([1e15,1e10,1])  # starting values (A, s, E) for the fit
# find optimal parameters 
params, cov=optimize.curve_fit(TLFOK,x_data,y_data,inis)
# params are the best fit values for the parameters 
# cov is the covariance of the best fit parameters 
plt.subplot(2,1, 1);  
plt.scatter(x_data, y_data, label='Glocanin TL #1');
plt.plot(x_data, TLFOK(x_data, *params),
c='r',linewidth=3, label='FOK-TL equation'); 
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('TL signal [a.u.]');
plt.xlabel(r'Temperature T [K]');
plt.xlim(375,550);
plt.subplot(2,1, 2); 
plt.plot(x_data,TLFOK(x_data, *params)-\
y_data,c='r',linewidth=2,label='Residuals');
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('Residuals');
plt.xlabel(r'Temperature T [K]');
plt.xlim(375,550);
plt.ylim(-0.0001,.0001);
plt.tight_layout()
A=format(params[0],"10.1E")
dA = format(np.sqrt(cov[0][0]),"10.1E")
s=format(params[1],"10.1E")
ds = format(np.sqrt(cov[1][1]),"10.1E")
E=round(params[2],3)
dE = round(np.sqrt(cov[2][2]),7)
res=TLFOK(x_data, *params)-y_data
FOM=round(100*np.sum(abs(res))/np.sum(y_data),3)
myTable = PrettyTable([ "A","s (s^-1)","ds","E(eV)","dE"]) 
myTable.add_row([A,s,ds,E,dE]);
print('FOM=',FOM, ' %')
print(myTable)
plt.show()