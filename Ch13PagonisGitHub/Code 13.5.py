#Fitting CW-IRSL signal with stretched exponential plus exponential
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
from prettytable import PrettyTable
from scipy.special import lambertw 
import warnings
data = np.loadtxt('t100.txt')
x_data,y_data=data[:,0], data[:,1]
y_data=y_data/max(y_data)
def STRETCH(t, A,tau,beta,tau2):
    CW=A*np.exp(-(t/tau)**beta)+ (1-A)*np.exp(-t/tau2) 
    return CW
def total_CW(t, *inis): 
    u=np.array([0 for i in range(len(x_data))])
    As, taus,betas,taus2=  inis
    for i in range(nPks):        
        u=u+STRETCH(t,As,taus,betas,taus2)
    return u
P=int(max(x_data)) 
t = np.linspace(0, P, P)
nPks=1
inis=[.5,10,.5,100]
params, cov = optimize.curve_fit(total_CW,\
x_data,y_data,p0=inis,maxfev=10000)   
plt.scatter(x_data, y_data,c='r',label='FL1 feldspar TR-IRSL ')
plt.plot(x_data, total_CW(x_data, 
 *params),c='black',label='Analytical equation')
for i in range(0,nPks): 
    CWi=STRETCH(t, *params)
    plt.plot(t,CWi)
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)  
plt.ylabel('CW-IRSL [a.u.]')
plt.xlabel('Time [$\mu$s]')
res=total_CW(x_data, *params)-y_data
FOM=100*np.sum(abs(res))/np.sum(y_data)
print('FOM=',round(FOM,1),' %')
As=round(params[0],2) 
As2=round(params[3],4) 
taus=round(params[1],1) 
taus2=int(params[3]) 
dAs=round(np.sqrt(cov[0][0]),4)
dtaus=round(np.sqrt(cov[1][1]),2)
dbeta=round(np.sqrt(cov[2][2]),3)
beta=round(params[2],3) 
dtaus2=int(np.sqrt(cov[3][3]))
myTable = PrettyTable([ "A (a.u.)",\
'tau (s)','dtau','beta','dbeta','tau2','dtau2']) 
myTable.add_row([As,taus,dtaus,beta,dbeta,taus2,dtaus2])
print(myTable)
plt.show()