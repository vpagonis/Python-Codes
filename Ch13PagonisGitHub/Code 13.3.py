# Analysis of TR-OSL experimental data in alumina
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
from prettytable import PrettyTable 
data = np.loadtxt('aluminaxy1.txt') 
x_data, y_data = np.array(data[:, 0][8:20]), \
np.array(data[:, 1][8:20])
kB=8.617e-5
def tq(x_data, N,c, W):
    return N /(1+c*np.exp(-W/(kB*(273+x_data))))
params, cov = optimize.curve_fit(tq,\
x_data, y_data,bounds=(0,[.04,1e12,1]))
plt.subplot(1,2, 1)
plt.plot(x_data,y_data,"o",c='black',label='Experimental')
plt.plot(x_data, tq(x_data, *params),
label='Best fit',c='r',linewidth=2)
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel(r'Luminescence Lifetime $\tau$ [s]')
plt.xlabel(r'Stimulation Temperature [$^{o}$C]')
plt.title('(a)   Finding C,W')
plt.text(90,.02,"Luminescence")
plt.text(90,.017,r"lifetime $\tau$")
C,W=format(params[1],"10.1E"),round(params[2],4)
dC = round(np.sqrt(cov[1][1]),4)
dW = round(np.sqrt(cov[2][2]),3)
data = np.loadtxt('aluminaxy2.txt') 
x0, y0 = np.array(data[:, 0]), np.array(data[:, 1])
x_data=1/(kB*(273+np.array(x0)[4:10]))
y_data=np.log(np.array(y0)[4:10])
def arrh(x_data,slope,inter):
    return slope*x_data+inter
params, cov = optimize.curve_fit(arrh,x_data, y_data)
plt.subplot(1,2, 2)
plt.plot(x_data,y_data,"o",c="black",label='Experimental')
x_vals=np.arange(30,35,.1)
plt.plot(x_vals, arrh(x_vals, *params),
label='Best fit',c='r',linewidth=2)
Eth=-round(params[0],3)
dEth = round(np.sqrt(cov[0][0]),3)
plt.title(r'(b)    Finding E$_{th}$')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel(r'ln($\tau$)')
plt.xlabel(r'1/(kT) [eV$^{-1}$]')
plt.text(31, 14.95,'Arrhenius plot')
plt.tight_layout()
myTable=PrettyTable(["C","dC","W (eV)","dW", "Eth (eV)","dEth"])  
myTable.add_row([C,dC,W,dW,Eth,dEth])
print(myTable)
plt.show()