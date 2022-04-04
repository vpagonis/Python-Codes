# LM-OSL deconvolution with KV-LM equation plus background
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
from prettytable import PrettyTable
from scipy.special import wrightomega 
import warnings
warnings.filterwarnings("ignore")
data = np.loadtxt('K120.txt')
x_data,y_data = data[:, 0], data[:, 1] 
y_data=y_data/max(y_data)
def KVLMOSL(t, A,c,sprime):
    zLMOSL=(1/c)-np.log(c)+sprime*t**2/(2*P)
    lam=wrightomega(zLMOSL)
    LMOSL=A*t/(lam+lam**2)   
    return LMOSL
def total_LMOSL(t, *inis): 
    u=np.array([0 for i in range(len(x_data))])
    As, cs ,sprimes=inis[0:nPks], inis[nPks:2*nPks],\
    inis[2*nPks:3*nPks]
    for i in range(nPks):        
        u=u+KVLMOSL(t,As[i],cs[i],sprimes[i])    
    ubgd=bgd*t/300
    u=u+ubgd
    return u
P=int(max(x_data)) 
nPks=1
t = np.linspace(0, P, P)
bgd=y_data[-1]-.07       # adjusted parameter for better fit
inis=[1,.1,1e-4]
params, cov = optimize.curve_fit(total_LMOSL,\
x_data,y_data,p0=inis,maxfev=10000)   
plt.scatter(x_data, y_data,c='r',label=r'LM-IRSL')
plt.plot(x_data, total_LMOSL(x_data, *params),c='black',\
label='KV-LM equation',linewidth=1)
plt.plot(t,bgd*t/P)
plt.plot(t,KVLMOSL(t,*params))
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)  
plt.ylabel('LM-IRSL [a.u.]')
plt.xlabel(r'Stimulation time [s]')
res=total_LMOSL(x_data, *params)-y_data
FOM=round(100*np.sum(abs(res))/np.sum(y_data),2)
print('FOM=',FOM, '%')
As=[round(x,3) for x in params[0:nPks]]
cs=[round(x,2) for x in params[nPks:2*nPks]]
sprimes=[round(x,4) for x in params[2*nPks:3*nPks]] 
dAs=[round(np.sqrt(cov[x][x]),2) for x in range(nPks)]
dcs=[round(np.sqrt(cov[x][x]),2) for x in range(nPks,2*nPks)]
dsprimes=[round(np.sqrt(cov[x][x]),5) for x in\
range(2*nPks,3*nPks)]
myTable = PrettyTable([ "A (a.u.)","dA",\
'c','dc','s (s^-1)','ds']) 
myTable.add_row([As[0],dAs[0],cs[0],dcs[0],sprimes[0],\
 dsprimes[0]])
print(myTable)
plt.show()