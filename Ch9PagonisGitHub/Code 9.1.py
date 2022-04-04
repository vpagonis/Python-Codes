#  deconvolution with FOK-CW equation KST4 feldspar CW-IRSL Data 
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
from prettytable import PrettyTable
import warnings
data = np.loadtxt('KST4ph300IR.txt')
x_data,y_data=data[:,0][1:800], data[:,1][1:800]
y_data=y_data/max(y_data)
def FOKCW(t, A,tau):
    CW=A*np.exp(-t/tau)   
    return CW
def total_CW(t, *inis): 
    u=np.array([0 for i in range(len(x_data))])
    As, taus=    inis[0:nPks], inis[nPks:2*nPks]
    bgd=inis[-1]
    for i in range(nPks):        
        u=u+FOKCW(t,As[i],taus[i])
    u=u+bgd 
    return u
P=int(max(x_data)) 
t = np.linspace(0, P, P)
nPks=2
inis=[.5,10,.5,50,.01]
params, cov = optimize.curve_fit(total_CW,\
x_data,y_data,p0=inis,maxfev=10000)   
plt.scatter(x_data, y_data,c='r',label='KST4 feldspar CW-IRSL ')
plt.plot(x_data, total_CW(x_data, 
 *params),c='black',label='FOK-CW equation',linewidth=1)
for i in range(0,nPks): 
    CWi=FOKCW(t, params[i],params[nPks+i])
    plt.plot(t,CWi)
plt.plot(t,[params[-1]]*len(t))
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)  
plt.ylabel('CW-IRSL [a.u.]')
plt.xlabel(r'Stimulation time [s]')
res=total_CW(x_data, *params)-y_data
FOM=100*np.sum(abs(res))/np.sum(y_data)
print('FOM=',round(FOM,1),' %')
As=[round(x,3) for x in params[0:nPks]]
taus=[round(x,1) for x in params[nPks:2*nPks]] 
dAs=[round(np.sqrt(cov[x][x]),3) for x in range(0,nPks)]
dtaus=[round(np.sqrt(cov[x][x]),2) for x in\
range(nPks,2*nPks)]
bgd=round(params[-1],3) 
myTable = PrettyTable([ "A (a.u.)","dA",\
'tau (s)','dtau (s)','bgd']) 
myTable.add_row([As[0],dAs[0],taus[0],dtaus[0],bgd])
myTable.add_row([As[1],dAs[1],taus[1],dtaus[1],' '])
print(myTable)
plt.show()