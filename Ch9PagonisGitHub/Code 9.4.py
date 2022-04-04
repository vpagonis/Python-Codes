# KST4 CW-IRSL deconvolution with GOK-CW equation 
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
from prettytable import PrettyTable
import warnings
warnings.filterwarnings("ignore")
data = np.loadtxt('KST4ph300IR.txt')
x_data,y_data=data[:,0][1:800], data[:,1][1:800]
y_data=y_data/max(y_data)
def GOKCW(t, A,b,sprime, bgd): 
    ITL=A*(1+sprime*(b-1)*t)**(-b/(b-1))  
    return ITL
def total_CW(t, *inis): 
    u=np.array([0 for i in range(len(x_data))])
    As, bs ,sprimes=inis[0:nPks], inis[nPks:2*nPks],\
    inis[2*nPks:3*nPks] 
    bgd=inis[-1] 
    for i in range(nPks):        
        u=u+GOKCW(t,As[i],bs[i],sprimes[i],bgd)    
    ubgd=bgd
    u=u+ubgd
    return u
P=int(max(x_data))
t = np.linspace(0, P,P)
nPks=1
inis=[1,1.5,.1,.01]
params, cov = optimize.curve_fit(total_CW,\
x_data,y_data,p0=inis,maxfev=10000)   
bgd=params[-1]
res=total_CW(x_data, *params)-y_data
FOM=round(100*np.sum(abs(res))/np.sum(y_data),2)
print('FOM=',FOM)
As=[round(x,2) for x in params[0:nPks]]
bs=[round(x,2) for x in params[nPks:2*nPks]]
sprime=[round(x,3) for x in params[2*nPks:3*nPks]] 
bgd=round(params[-1],3) 
dAs=[round(np.sqrt(cov[x][x]),3) for x in range(nPks)]
dbs=[round(np.sqrt(cov[x][x]),3) for x in range(nPks,2*nPks)]
dsprime=[round(np.sqrt(cov[x][x]),4)\
for x in range(2*nPks,3*nPks)]
myTable = PrettyTable([ "A (a.u.)","dA (a.u)",\
'b','db',"s' (s^-1)","ds'"]) 
myTable.add_row([As[0],dAs[0],bs[0],dbs[0],sprime[0],\
dsprime[0]])
print(myTable)