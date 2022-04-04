# Isothermal analysis for Durango apatite
#  deconvolution of 220degC ITL data with sum of 3 exponentials 
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
from prettytable import PrettyTable
import warnings
warnings.filterwarnings("ignore")
data = np.loadtxt('DurangoITL220.txt')
x_data,y_data = data[:, 0][3:1000], data[:, 1][3:1000]
y_data=y_data/max(y_data)
def oneexpon(x_data,N,tau):
    u=np.abs(N)*(np.exp(-x_data/np.abs(tau)))
    return u
def total_ITL(t, *inis): 
    u=np.array([0 for i in range(len(x_data))])
    Ns, taus =    inis[0:nPks], inis[nPks:2*nPks]
    for i in range(nPks):        
        u=u+oneexpon(t,Ns[i],taus[i])
    return u
t = np.linspace(0, len(x_data), len(x_data))
nPks=3
N=[1,.5,.1]
tau=[10,100,1000]
inis=N+tau
params, cov = optimize.curve_fit(total_ITL,\
x_data,y_data,p0=inis,maxfev=10000)   
plt.scatter(x_data, y_data,c='r',label='Durango ITL signal')
plt.plot(x_data, total_ITL(x_data, 
 *params),c='black',label='Sum of exponentials',linewidth=1)
sums=[0]*nPks
for i in range(0,nPks): 
    ITLi=oneexpon(t, params[i],params[nPks+i])
    sums[i]=np.sum(ITLi)
    plt.plot(t,ITLi)
pc1=round(100*sums[0]/(sums[0]+sums[1]+sums[2]),1) 
pc2=round(100*sums[1]/(sums[0]+sums[1]+sums[2]),1) 
pc3=round(100*sums[2]/(sums[0]+sums[1]+sums[2]),1)
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)  
plt.ylabel('ITL [a.u.]')
plt.xlabel(r'Stimulation time [s]')
res=total_ITL(x_data, *params)-y_data
FOM=100*np.sum(abs(res))/np.sum(y_data)
Ns=[round(x,3) for x in params[0:nPks]]
taus=[round(x,1) for x in params[nPks:2*nPks]] 
dN=[round(np.sqrt(cov[x][x]),3) for x in range(3)]
dtaus=[round(np.sqrt(cov[x][x]),1) for x in range(3,6)]
myTable = PrettyTable([ "Curve","N (a.u.)","dN (a.u)",\
'tau (s)',"dtau (s)","%"]) 
myTable.add_row(["1",Ns[0],dN[0],taus[0],dtaus[0],pc1])
myTable.add_row(["2",Ns[1],dN[1],taus[1],dtaus[1],pc2])
myTable.add_row(["3",Ns[2],dN[2],taus[2],dtaus[2],pc3])
print('FOM=',round(FOM,2),' %')
print(myTable)
plt.show()