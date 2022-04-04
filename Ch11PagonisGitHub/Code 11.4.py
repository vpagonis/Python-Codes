#Fitting CW-IRSL signal with stretched exponential
#  deconvolution of feldspar CW-IRSL with stretched exponential 
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
from prettytable import PrettyTable
import warnings
data = np.loadtxt('KST4ph300IR.txt')
x_data,y_data=data[:,0][0:100], data[:,1][0:100]
y_data=y_data/max(y_data)
def STRETCH(t, A,tau,beta):
    CW=A*np.exp(-(t/tau)**beta)   
    return CW
def total_CW(t, *inis): 
    u=np.array([0 for i in range(len(x_data))])
    As, taus=    inis[0:nPks], inis[nPks:2*nPks]
    beta=inis[-1]
    for i in range(nPks):        
        u=u+STRETCH(t,As[i],taus[i],beta)
    return u
P=int(max(x_data)) 
t = np.linspace(0, P, P)
nPks=1
inis=[1,5,.1]
params, cov = optimize.curve_fit(total_CW,\
x_data,y_data,p0=inis,maxfev=10000)   
plt.scatter(x_data, y_data,c='r',label='KST4 feldsparCW-IRSL ')
plt.plot(x_data, total_CW(x_data, 
 *params),c='black',label='Stretched Exponential equation')
for i in range(0,nPks): 
    CWi=STRETCH(t, params[i],params[nPks+i],params[-1])
    plt.plot(t,CWi)
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)  
plt.ylabel('CW-IRSL [a.u.]')
plt.xlabel(r'Stimulation time [s]')
res=total_CW(x_data, *params)-y_data
FOM=100*np.sum(abs(res))/np.sum(y_data)
As=[round(x,2) for x in params[0:nPks]]
taus=[round(x,1) for x in params[nPks:2*nPks]] 
dAs=[round(np.sqrt(cov[x][x]),2) for x in range(0,nPks)]
dtaus=[round(np.sqrt(cov[x][x]),2) for x in\
range(nPks,2*nPks)]
dbeta=round(np.sqrt(cov[2][2]),2)
beta=round(params[-1],2) 
myTable = PrettyTable([ "A (a.u.)","dA",\
'tau (s)','dtau (s)','beta','dbeta']) 
myTable.add_row([As[0],dAs[0],taus[0],dtaus[0],beta,dbeta])
print('FOM=',round(FOM,1),' %')
print(myTable)
plt.show()