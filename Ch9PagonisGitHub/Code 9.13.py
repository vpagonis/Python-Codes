# Deconvolution with the KP-LM equation, fixed bgd
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
from prettytable import PrettyTable 
import warnings
warnings.filterwarnings("ignore")
data = np.loadtxt('K120.txt')
x_data,y_data=data[:,0], data[:,1]
y_data=y_data/max(y_data)
def KVLMOSL(x, A,sprime,rho):
    F=np.log(1+sprime*x**2/(2*P))
    LMOSL= abs(A)*np.exp (-rho*\
	F** 3.0)*(F**2.0)/(1+sprime*x**2/(2*P))  
    return LMOSL
def total_LMOSL(x, *inis): 
    u=np.array([0 for i in range(len(x_data))])
    As, sprimes=inis[0:nPks], inis[nPks:2*nPks]
    rho=inis[-1]
    for i in range(nPks):        
        u=u+KVLMOSL(x,As[i],sprimes[i],rho)    
    return u
nPks=2
P=int(max(x_data))
inis=[.5,10,.5,100,.005]
params, cov = optimize.curve_fit(total_LMOSL,\
x_data, y_data,p0=inis)
plt.scatter(x_data, y_data, label='Experiment')
plt.plot(x_data, total_LMOSL(x_data, *params),
         label='KP-LM equation',linewidth=3)
for i in range(0,nPks): 
    FOKLMi=KVLMOSL(x_data, params[i],params[nPks+i],params[-1])
    plt.plot(x_data,FOKLMi)
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('LM-IRSL signal [a.u.]')
plt.xlabel('Time [s]')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
A1,sprime1,A2,sprime2,rho=[round(params[x],3) for x in range(5)]
dA1,dsprime1,dA2,dsprime2,drho=[round(np.sqrt(cov[x][x]),4)\
for x in range(5)]
res=total_LMOSL(x_data, *params)-y_data
FOM=round(100*np.sum(abs(res))/np.sum(y_data),2)
myTable=PrettyTable(["A",'dA', "rho",  "d(rho)",\
"s'(s^-1)","ds","FOM"])  
myTable.add_row([A1,dA1,rho,drho, sprime1, dsprime1,FOM])
myTable.add_row([A2,dA2,' ',' ', sprime2, dsprime2,' '])
print(myTable) 
plt.show()