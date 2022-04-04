#Analysis of 2-component LM-OSL signal usig FOK-LM eqt
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
from prettytable import PrettyTable
import warnings
warnings.filterwarnings("ignore")
data = np.loadtxt('K120.txt')
x_data,y_data = data[:, 0], data[:, 1] 
y_data =y_data/max(y_data)
def LM(x_data,N,xmax):
    u=1.6487*np.abs(N)*(x_data/abs(xmax))*(np.exp(-(x_data**\
    2.0)/(2*(abs(xmax)**2.0))))
    return u
def total_FOKLM(t, *inis): 
    u=np.array([0 for i in range(len(x_data))])
    Ns, xmaxs =    inis[0:nPks], inis[nPks:2*nPks]
    for i in range(nPks):        
        u=u+LM(t,Ns[i],xmaxs[i])
    u=u+bgd*t/P
    return u
nPks= 2 
P=int(max(x_data))
t=np.linspace(0,P,P)
inis=[50,3,100,1]
bgd=y_data[-1]
params,cov =optimize.curve_fit(total_FOKLM,x_data,\
y_data,p0=inis)
params, cov = optimize.curve_fit(total_FOKLM,\
x_data,y_data,p0=inis,maxfev=10000)   
plt.scatter(x_data, y_data,c='r',label=r'LM-IRSL')
plt.plot(x_data, total_FOKLM(x_data, 
 *params),c='black',label='Transformed FOK-LM eqt',linewidth=1)
plt.plot(t,bgd*t/P)
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)  
plt.ylabel('LM-IRSL [a.u.]')
plt.xlabel(r'Stimulation time [s]')
for i in range(0,nPks): 
    FOKLMi=LM(x_data, params[i],params[nPks+i])
    plt.plot(x_data,FOKLMi)
res=total_FOKLM(x_data, *params)-y_data
FOM=round(100*np.sum(abs(res))/np.sum(y_data),1)
print('FOM=',FOM,' %')
plt.show()
Ns=[round(x,2) for x in params[0:nPks]]
xmaxs=[round(abs(x),2) for x in params[nPks:2*nPks]] 
dN=[round(np.sqrt(cov[x][x]),2) for x in range(nPks)]
dxmaxs=[round(np.sqrt(cov[x][x]),2) for x in  range(nPks,2*nPks)]
myTable = PrettyTable([ "Im (a.u.)","dIm (a.u)",\
'tm (s)',"dtm (s)"]) 
for j in range(nPks):
    myTable.add_row([Ns[j],dN[j],xmaxs[j],dxmaxs[j]])
print(myTable)
