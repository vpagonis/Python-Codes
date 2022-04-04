#  deconvolution with KV-ITL equation Durango 220deg ITL Data 
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
from prettytable import PrettyTable
from scipy.special import wrightomega 
import warnings
data = np.loadtxt('DurangoITL220.txt')
x_data,y_data = data[:, 0][3:1000], data[:, 1][3:1000]
y_data=y_data/max(y_data)
plt.plot(x_data,y_data)
def KVITL(t, A,c,sprime):
    zITL=(1/c)-np.log(c)+sprime*t
    lam=wrightomega(zITL)
    ITL=A/(lam+lam**2)   
    return ITL
def total_ITL(t, *inis): 
    u=np.array([0 for i in range(len(x_data))])
    As, cs ,sprimes=    inis[0:nPks], inis[nPks:2*nPks],\
    inis[2*nPks:3*nPks]
    bgd=inis[-1]
    for i in range(nPks):        
        u=u+KVITL(t,As[i],cs[i],sprimes[i])
    u=u+bgd 
    return u
t = np.linspace(1, 1000, 1000)
nPks=1
A=[max(y_data)]*nPks
lowA, highA=[0.01*x for x in A], [200 for x in A]
c=[1]*nPks
lowc, highc= [0.001*x for x in c], [1e4*x for x in c]
sprime=[.01]*nPks 
lowsprime, highsprime=  [0.001*x for x in sprime],\
[10*x for x in sprime]
bgd, lowbgd, highbgd=[.1,0,.15]
inis=A+c+sprime+[bgd]
lowbnds=lowA+lowc+lowsprime+[lowbgd]
highbnds=highA+highc+highsprime+[highbgd]
params, cov = optimize.curve_fit(total_ITL,\
x_data,y_data,p0=inis,bounds=(lowbnds,highbnds),maxfev=10000)   
plt.scatter(x_data, y_data,c='r',label='Durango ITL signal')
plt.plot(x_data, total_ITL(x_data, 
 *params),c='black',label='KV-ITL',linewidth=1)
sums=[0]*nPks
for i in range(0,nPks): 
    ITLi=KVITL(t, params[i],params[nPks+i], params[2*nPks+i])
    sums[i]=np.sum(ITLi)
    plt.plot(t,ITLi)
bgdarea=max(t)*params[-1]
plt.plot(t,[params[-1]]*len(t))
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)  
plt.ylabel('ITL-OSL [a.u.]')
plt.xlabel(r'Stimulation time [s]')
res=total_ITL(x_data, *params)-y_data
FOM=100*np.sum(abs(res))/np.sum(y_data)
print('FOM=',round(FOM,1),' %')
As=[round(x,1) for x in params[0:nPks]]
cs=[round(x,2) for x in params[nPks:2*nPks]]
sprimes=[round(x,3) for x in params[2*nPks:3*nPks]] 
dAs=[round(np.sqrt(cov[x][x]),1) for x in range(nPks)]
dcs=[round(np.sqrt(cov[x][x]),2) for x in range(nPks,2*nPks)]
dsprimes=[round(np.sqrt(cov[x][x]),3) for x in\
range(2*nPks,3*nPks)]
bgd=round(params[-1],3) 
myTable = PrettyTable([ "A (a.u.)","dA",\
'c','dc',"s' (s^-1)","ds'",'bgd']) 
myTable.add_row([As[0],dAs[0],cs[0],dcs[0],\
sprimes[0],dsprimes[0],bgd])
print(myTable)
plt.show()