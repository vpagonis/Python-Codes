#  deconvolution with KV-CW equation KST4 feldspar CW-IRSL Data 
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
from prettytable import PrettyTable
from scipy.special import wrightomega 
import warnings
data = np.loadtxt('KST4ph300IR.txt')
x_data,y_data=data[:,0][1:800], data[:,1][1:800]
y_data=y_data/max(y_data)
def KVCW(t, A,c,sprime):
    zCW=(1/c)-np.log(c)+sprime*t
    lam=wrightomega(zCW)
    CW=A/(lam+lam**2)   
    return CW
def total_CW(t, *inis): 
    u=np.array([0 for i in range(len(x_data))])
    As, cs ,sprimes=    inis[0:nPks], inis[nPks:2*nPks],\
    inis[2*nPks:3*nPks]
    bgd=inis[-1]
    for i in range(nPks):        
        u=u+KVCW(t,As[i],cs[i],sprimes[i])
    u=u+bgd 
    return u
t = np.linspace(0, 200, 200)
nPks=1
A=[max(y_data)]*nPks
lowA, highA=[0.01*x for x in A], [200 for x in A]
c=[1]*nPks
lowc, highc= [0.001*x for x in c], [1e4*x for x in c]
sprime=[.1]*nPks 
lowsprime, highsprime=  [0.001*x for x in sprime],\
[100*x for x in sprime]
bgd, lowbgd, highbgd=[.1,0,.15]
inis=A+c+sprime+[bgd]
lowbnds=lowA+lowc+lowsprime+[lowbgd]
highbnds=highA+highc+highsprime+[highbgd]
params, cov = optimize.curve_fit(total_CW,\
x_data,y_data,p0=inis,bounds=(lowbnds,highbnds),maxfev=10000)   
plt.scatter(x_data, y_data,c='r',label='KST4 feldspar CW-IRSL ')
plt.plot(x_data, total_CW(x_data, 
 *params),c='black',label='KV-CW equation',linewidth=1)
for i in range(0,nPks): 
    CWi=KVCW(t, params[i],params[nPks+i], params[2*nPks+i])
    plt.plot(t,CWi)
plt.plot(t,[params[-1]]*len(t))
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)  
plt.ylabel('CW-IRSL [a.u.]')
plt.xlabel(r'Stimulation time [s]')
res=total_CW(x_data, *params)-y_data
FOM=100*np.sum(abs(res))/np.sum(y_data)
print('FOM=',round(FOM,1),' %')
As=[round(x,2) for x in params[0:nPks]]
cs=[round(x,2) for x in params[nPks:2*nPks]]
sprimes=[round(x,2) for x in params[2*nPks:3*nPks]] 
dAs=[round(np.sqrt(cov[x][x]),2) for x in range(nPks)]
dcs=[round(np.sqrt(cov[x][x]),2) for x in range(nPks,2*nPks)]
dsprimes=[round(np.sqrt(cov[x][x]),2) for x in\
range(2*nPks,3*nPks)]
bgd=round(params[-1],2) 
myTable = PrettyTable([ "A (a.u.)","dA",\
'c','dc','s (s^-1)','ds (s^-1)']) 
myTable.add_row([As[0],dAs[0],cs[0],dcs[0],sprimes[0],\
dsprimes[0]])
print(myTable)
plt.show()