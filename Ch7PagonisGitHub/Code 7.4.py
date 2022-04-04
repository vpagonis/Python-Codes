#  deconvolution with GOK-ITL equation Durango 220deg ITL Data 
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
from prettytable import PrettyTable
import warnings
warnings.filterwarnings("ignore")
data = np.loadtxt('DurangoITL220.txt')
x_data,y_data = data[:, 0][3:1000], data[:, 1][3:1000]
y_data=y_data/max(y_data)
def GOKITL(t, A, b, sprime):
    ITL=A*(1+(b-1)*sprime*t)**(-b/(b-1))
    return ITL
def total_GOKITL(t, *inis): 
    u=np.array([0 for i in range(len(x_data))])
    As, bs ,sprimes=inis[0:nPks], inis[nPks:2*nPks],\
    inis[2*nPks:3*nPks] 
    bgd=inis[-1] 
    for i in range(nPks):        
        u=u+GOKITL(t,As[i],bs[i],sprimes[i])          
    ubgd=bgd
    u=u+ubgd
    return u
t = np.linspace(0, 1000, 1000)
nPks=1
A=[max(y_data)]*nPks
lowA, highA=[0.01*x for x in A], [1 for x in A]
b=[1.5]*nPks
lowb, highb= [1.001 for x in b], [2 for x in b]
sprime=[10]*nPks 
lowsprime, highsprime=  [1e-4 for x in sprime],\
[1e4 for x in sprime]
bgd=.01
lowbgd,highbgd=0,.3
inis=A+b+sprime+[bgd]
lowbnds=lowA+lowb+lowsprime+[lowbgd]
highbnds=highA+highb+highsprime+[highbgd]
params, cov = optimize.curve_fit(total_GOKITL,\
x_data,y_data,p0=inis,bounds=(lowbnds,highbnds),maxfev=10000)   
res=total_GOKITL(x_data, *params)-y_data
FOM=100*np.sum(abs(res))/np.sum(y_data)
As=[round(x,3) for x in params[0:nPks]]
bs=[round(x,1) for x in params[nPks:2*nPks]]
sprimes=[round(x,4) for x in params[2*nPks:3*nPks]] 
dAs=[round(np.sqrt(cov[x][x]),3) for x in range(nPks)]
dbs=[round(np.sqrt(cov[x][x]),1) for x in range(nPks,2*nPks)]
dsprimes=[round(np.sqrt(cov[x][x]),4) for x in\
range(2*nPks,3*nPks)]
bgd=round(params[-1],3) 
print('FOM=',round(FOM,1)," %")
myTable = PrettyTable(["A (a.u.)","dA",\
'b','db',"s' (s^-1)","ds' (s^-1)"]) 
myTable.add_row([As[0],dAs[0],bs[0],dbs[0],\
sprimes[0],dsprimes[0]])
print(myTable)
