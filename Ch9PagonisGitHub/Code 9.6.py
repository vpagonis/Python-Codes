# Analysis of 3-component LM-OSL signal
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
from prettytable import PrettyTable
import warnings
warnings.filterwarnings("ignore")
data = np.loadtxt('CaF2LMOSL.txt')
x_data,y_data = data[:, 0], data[:, 1] 
def LM(x,N,tau):
    u=np.abs(N)*(x/P)*(np.exp(-(x**2.0)\
    /(2*P*abs(tau))))
    return u
def total_FOKLM(x, *inis): 
    u=np.array([0 for i in range(len(x_data))])
    Ns, taus =    inis[0:nPks], inis[nPks:2*nPks]
    for i in range(nPks):        
        u=u+LM(x,Ns[i],taus[i])
    u=u+bgd*x/P
    return u
nPks= 3 
P=int(max(x_data))
t=np.linspace(0,P,P)
inis=[1400,1,800,.1,500,.01]
bgd=y_data[-1]
params,cov =optimize.curve_fit(total_FOKLM,x_data,\
y_data,p0=inis)
params, cov = optimize.curve_fit(total_FOKLM,\
x_data,y_data,p0=inis,maxfev=10000)   
plt.scatter(x_data, y_data,c='r',label=r'CaF$_2$:Mn  LM-OSL')
plt.plot(x_data, total_FOKLM(x_data, 
 *params),c='black',label='Original FOK-LM eqt',linewidth=1)
totalArea=sum(total_FOKLM(x_data,  *params))
sums,pc=[0]*nPks, [0]*nPks
for i in range(0,nPks): 
    FOKLMi=LM(x_data, params[i],params[nPks+i])
    sums[i]=np.sum(FOKLMi)
    plt.plot(x_data,FOKLMi)
plt.plot(t,bgd*t/P)
for j in range(nPks):
    pc[j]=round(100*sums[j]/totalArea,1) 
pcbgd=round(100*sum(bgd*x_data/P)/totalArea)
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)  
plt.ylabel('LM-OSL [a.u.]')
plt.xlabel(r'Stimulation time [s]')
res=total_FOKLM(x_data, *params)-y_data
FOM=round(100*np.sum(abs(res))/np.sum(y_data),1)
print('FOM=',FOM,' %')
plt.show()
Ns=[round(x,1) for x in params[0:nPks]]
taus=[round(x,2) for x in params[nPks:2*nPks]] 
dN=[round(np.sqrt(cov[x][x]),1) for x in range(3)]
dtaus=[round(np.sqrt(cov[x][x]),2) for x in range(3,6)]
myTable = PrettyTable([ "N (a.u.)","dN (a.u)",\
'tau (s)',"dtau (s)","Area [%]"]) 
for j in range(nPks):
    myTable.add_row([Ns[j],dN[j],taus[j],dtaus[j],pc[j]])
myTable.add_row(['','','','','bgd='+str(pcbgd)+'%'])
print(myTable)