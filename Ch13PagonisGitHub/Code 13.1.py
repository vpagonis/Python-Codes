# Analysis of TR-OSL in sedimentary quartz
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
from scipy import optimize
from prettytable import PrettyTable 
data = np.loadtxt('chithamboTROSLqzdata.txt')
x_data,y_data = data[:, 0][0:27], data[:, 1][0:27]
plt.plot(x_data,y_data,"o")
def SEfit(x_data,N,Dc):
    u=N*(1-np.exp(-x_data/np.abs(Dc)))
    return u
init_vals=[200,40]
params, cov = optimize.curve_fit(SEfit,\
x_data, y_data,p0=init_vals)
dN = round(np.sqrt(cov[0][0]),1)
dDo = round(np.sqrt(cov[1][1]),1)
x_vals=np.arange(0,60,1)
plt.plot(x_vals, SEfit(x_vals, *params[0:2]),c="b",\
label='SE fit')     
plt.scatter(x_data, y_data, label='Experiment')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylim(0,260)
plt.ylabel('TR-OSL [a.u.]')
plt.xlabel('Time [$\mu$s]')
plt.title('TR-OSL sedimentary quartz')
plt.text(200,150,'470 nm LEDs')
plt.tight_layout()
myTable = PrettyTable(["N", "dN","tau (microsec.)",\
"d(tau) (microsec.)"]) 
myTable.add_row([round(params[0],1),dN,round(params[1],1),dDo])
x_data,y_data = data[:, 0][28:99], data[:, 1][28:99]
plt.plot(x_data,y_data,"o")
def decayfit(x_data,N,Dc):
    u=N*(np.exp(-x_data/np.abs(Dc)))
    return u
init_vals=[200,40]
params, cov = optimize.curve_fit(decayfit,\
x_data, y_data,p0=init_vals)
dN = round(np.sqrt(cov[0][0]),1)
dDo = round(np.sqrt(cov[1][1]),1)
x_vals=np.arange(60,300,1)
plt.plot(x_vals, decayfit(x_vals, *params[0:2]),c="b")
plt.scatter(x_data, y_data)
plt.tight_layout()
myTable.add_row([round(params[0],1),dN,round(params[1],1),dDo])
print(myTable)
plt.show()