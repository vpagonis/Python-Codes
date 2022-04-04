#  deconvolution with MOK-ITL equation Durango 220deg ITL Data 
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
from prettytable import PrettyTable
import warnings
warnings.filterwarnings("ignore")
data = np.loadtxt('DurangoITL220.txt')
x_data,y_data = data[:, 0][3:1000], data[:, 1][3:1000]
y_data=y_data/max(y_data)
# b in this code represents the MOK parameter alpha
def MOKITL(t, A, b, sprime):
    F=np.exp(sprime*t)
    ITL=A*F/((F-b)**2)  
    return ITL
def total_MOKITL(t, *inis): 
    u=np.array([0 for i in range(len(x_data))])
    As, bs ,sprimes=inis[0:nPks], inis[nPks:2*nPks],\
    inis[2*nPks:3*nPks] 
    bgd=inis[-1] 
    for i in range(nPks):        
        u=u+MOKITL(t,As[i],bs[i],sprimes[i])          
    ubgd=bgd
    u=u+ubgd
    return u
t = np.linspace(0, 1000, 1000)
nPks=1
A=[100]*nPks
lowA, highA=[0.1 for x in A], [1e4 for x in A]
b=[0.1]*nPks
lowb, highb= [0.001 for x in b], [1 for x in b]
sprime=[.01]*nPks 
lowsprime, highsprime=  [1e-4 for x in sprime],\
[1e4 for x in sprime]
bgd=.01
lowbgd,highbgd=0,.3
inis=A+b+sprime+[bgd]
lowbnds=lowA+lowb+lowsprime+[lowbgd]
highbnds=highA+highb+highsprime+[highbgd]
params, cov = optimize.curve_fit(total_MOKITL,\
x_data,y_data,p0=inis,bounds=(lowbnds,highbnds),maxfev=10000)   
res=total_MOKITL(x_data, *params)-y_data
FOM=100*np.sum(abs(res))/np.sum(y_data)
print('FOM=',round(FOM,1),' %')
As=[round(x,3) for x in params[0:nPks]]
alpha=[round(x,2) for x in params[nPks:2*nPks]]
sprimes=[round(x,4) for x in params[2*nPks:3*nPks]] 
dAs=[round(np.sqrt(cov[x][x]),3) for x in range(nPks)]
dalpha=[round(np.sqrt(cov[x][x]),2) for x in range(nPks,2*nPks)]
dsprimes=[round(np.sqrt(cov[x][x]),4) for x in\
range(2*nPks,3*nPks)]
bgd=round(params[-1],2) 
myTable = PrettyTable(["A (a.u.)","dA",\
'alpha','dalpha',"s' (s^-1)","ds' (s^-1)"]) 
myTable.add_row([As[0],dAs[0],alpha[0],dalpha[0],\
sprimes[0],dsprimes[0]])
print(myTable)