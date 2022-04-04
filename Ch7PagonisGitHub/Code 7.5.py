# deconvolution with KP-ITL equation Durango 220deg ITL Data 
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
from prettytable import PrettyTable
import warnings
warnings.filterwarnings("ignore")
data = np.loadtxt('DurangoITL220.txt')
x_data,y_data = data[:, 0][3:1000], data[:, 1][3:1000]
y_data=y_data/max(y_data)
plt.subplot(2,1, 1)
def test_func(x, imax_fit,rho_fit, A_fit):
    return imax_fit*np.exp (-rho_fit*(np.log(1 + A_fit*x))\
    ** 3.0)*(np.log(1+A_fit*x)**2.0)/(1+x*A_fit)
params, cov = optimize.curve_fit(test_func,\
x_data, y_data)
drho= round(np.sqrt(cov[1][1]),5)
dA = round(np.sqrt(cov[2][2]),3)
plt.plot(np.log(x_data), y_data,'+',c='r', label='Experiment')
plt.plot(np.log(x_data), test_func(x_data, *params[0:4]),
label='KP-ITL equation')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('ITL signal [a.u.]')
plt.xlabel('ln(time)')
plt.text(1,.5,'Durango ITL signal')
plt.subplot(2,1, 2)
plt.plot(np.log(x_data),test_func(x_data, *params[0:4])-\
y_data,"o",label='Residuals')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('Residuals')
plt.xlabel('ln(time)')
plt.ylim(-.1,.1)
plt.tight_layout()
imax,rho, A=int(params[0]),round(params[1],5),\
round(params[2],2)
res=test_func(x_data, *params)-y_data
FOM=round(100*np.sum(abs(res))/np.sum(y_data),2)
myTable=PrettyTable(["imax", "rho",  "d(rho)",\
"sprime(s^-1)","dsprime","FOM"])  
myTable.add_row([imax,rho,drho, A, dA,FOM])
print(myTable) 
plt.show()