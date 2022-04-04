# Deconvolution with the KP-LM equation, fixed bgd
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
from prettytable import PrettyTable 
import warnings
warnings.filterwarnings("ignore")
data = np.loadtxt('J1000 50degC LMIRSL.txt')
x_data,y_data=data[:,0], data[:,1]
y_data=y_data/max(y_data)
def KVLMOSL(x, A,sprime,rho):
    F=np.log(1+abs(sprime)*x**2/(2*P))
    LMOSL= abs(A)*np.exp (-abs(rho)*\
	F** 3.0)*(F**2.0)/(1+abs(sprime)*x**2/(2*P))  
    return LMOSL
def total_LMOSL(x, *inis): 
    u=np.array([0 for i in range(len(x_data))])
    As, sprimes=inis[0:nPks], inis[nPks:2*nPks]
    rho=inis[-1]
    for i in range(nPks):        
        u=u+KVLMOSL(x,As[i],sprimes[i],rho)        
    u=u+bgd*x/P
    return u
nPks=2
P=int(max(x_data))
bgd=y_data[-1]
inis=[1,1,3,100,.01]
params, cov = optimize.curve_fit(total_LMOSL,\
x_data, y_data,p0=inis)
plt.scatter(x_data, y_data, label='J1000 feldspar')
plt.plot(x_data, total_LMOSL(x_data, *params),'r-',
         label='KP-LM equation',linewidth=3)
for i in range(0,nPks): 
    FOKLMi=KVLMOSL(x_data, params[i], params[i+nPks],
    params[-1])
    plt.plot(x_data,FOKLMi)
plt.plot(x_data,x_data*bgd/P)
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('LM-IRSL signal [a.u.]')
plt.xlabel('Time [s]')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
A1,sprime1,A2,sprime2=[round(params[x],2) for x in range(4)]
rho, drho=f'{params[-1]:.4f}',f'{np.sqrt(cov[4][4]):.4f}'
dA1,dsprime1,dA2,dsprime2=[round(np.sqrt(cov[x][x]),2)\
for x in range(4)]
res=total_LMOSL(x_data, *params)-y_data
FOM=round(100*np.sum(abs(res))/np.sum(y_data),2)
myTable=PrettyTable(["A",'dA', "rho",  "d(rho)",\
"s'(s^-1)","ds","FOM"])  
myTable.add_row([A1,dA1,rho,drho, sprime1, dsprime1,FOM])
myTable.add_row([A2,dA2,' ',' ', sprime2, dsprime2,' '])
print(myTable) 
plt.show()