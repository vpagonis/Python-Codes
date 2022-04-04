# Analysis of TR-OSL in microcline feldspar quartz
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
from scipy import optimize
from prettytable import PrettyTable 
data = np.loadtxt('FL1ONdata.txt') 
z=1.8
x_data, y_data = np.array(data[:, 0]), np.array(data[:, 1])
def TRfit(x_data, imax, rho, s):
    u=imax*(1-np.exp (-rho*(np.log(1 +z* s* x_data)) **3.0))
    return u
init_vals=[1,.007,10]
params, cov = optimize.curve_fit(TRfit, x_data, y_data, \
p0=init_vals)
plt.subplot(1,2, 1)
x_vals=np.arange(0,50,1)
plt.plot(x_vals, TRfit(x_vals, *params),c="black",\
label='Analytical fit')       
plt.scatter(x_data, y_data, c="g",label='Experiment')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylim(0,1.2)
plt.ylabel('TR-IRSL [a.u.]')
plt.xlabel('Time [$\mu$s]')
plt.title('(a) LED ON period')
imax=round(params[0],1)
rho=round(params[1],3)
s=params[2]
x_vals=np.arange(0,105,1)
data = np.loadtxt('FL1OFFdata.txt') 
x_data, y_data = np.array(data[:, 0]), np.array(data[:, 1])
u=imax*((np.exp (-rho*(np.log(1 +z* s* x_vals)) **3.0))-\
(np.exp (-rho*(np.log(1 +z* s* (x_vals+50)) **3.0))))+.02
plt.subplot(1,2, 2)
plt.plot(x_vals, u,c="green",\
label='Analytical fit')       
plt.scatter(x_data, y_data, c='r',label='Experiment')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylim(0,1.2)
plt.ylabel('TR-IRSL [a.u.]')
plt.xlabel('Time [$\mu$s]')
plt.title('(b) LED OFF period')
plt.text(20,0.5,'TR-IRSL microcline')
plt.tight_layout()
s=format(params[2]*1e6,'10.1E')
ds = format(np.sqrt(cov[2][2])*1e6,'10.1E')
drho = round(np.sqrt(cov[1][1]),3)
myTable=PrettyTable(["LED","imax", "rho",  "d(rho)",\
"s (s^-1)","ds (s^-1)"])  
myTable.add_row(["ON",imax,rho,drho, s, ds])
print(myTable)
plt.show()