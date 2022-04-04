# KST4 CW-IRSL deconvolution with MOK-CW equation 
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
from prettytable import PrettyTable
import warnings
warnings.filterwarnings("ignore")
data = np.loadtxt('KST4ph300IR.txt')
x_data,y_data=data[:,0][1:800], data[:,1][1:800]
y_data=y_data/max(y_data)
def MOKCW(t, A,b,sprime, bgd): 
    F=np.exp(sprime*t)
    ITL=A*F/((F-b)**2)  
    return ITL
def total_CW(t, *inis): 
    u=np.array([0 for i in range(len(x_data))])
    As, bs ,sprimes=inis[0:nPks], inis[nPks:2*nPks],\
    inis[2*nPks:3*nPks],
    bgd=inis[-1] 
    for i in range(nPks):        
        u=u+MOKCW(t,As[i],bs[i],sprimes[i],bgd)    
    ubgd=bgd
    u=u+ubgd
    return u
P=int(max(x_data)) 
t = np.linspace(0, P, P)
nPks=1
A=[max(y_data)/2]*nPks
lowA, highA=[0.01*x for x in A], [2*x for x in A]
b=[0.1]*nPks
lowb, highb= [0.001 for x in b], [1 for x in b]
sprime=[.01]*nPks 
lowsprime, highsprime=  [1e-4 for x in sprime], \
[1e4 for x in sprime]
bgd=.01
lowbgd,highbgd=0,.3
inis=A+b+sprime+[bgd]
lowbnds=lowA+lowb+lowsprime+[lowbgd]
highbnds=highA+highb+highsprime+[highbgd]
params, cov = optimize.curve_fit(total_CW,\
x_data,y_data,p0=inis,bounds=(lowbnds,highbnds),maxfev=10000)   
res=total_CW(x_data, *params)-y_data
FOM=round(100*np.sum(abs(res))/np.sum(y_data),2)
bgd=params[-1]
As=[round(x,3) for x in params[0:nPks]]
cs=[round(x,3) for x in params[nPks:2*nPks]]
sprime=[round(x,3) for x in params[2*nPks:3*nPks]] 
dAs=[round(np.sqrt(cov[x][x]),3) for x in range(nPks)]
dcs=[round(np.sqrt(cov[x][x]),3) for x in range(nPks,2*nPks)]
dsprime=[round(np.sqrt(cov[x][x]),3) \
for x in range(2*nPks,3*nPks)]
myTable = PrettyTable([ "A (a.u.)","dA (a.u)",\
'alpha','dalpha',"s' (s^-1)","ds'"]) 
myTable.add_row([As[0],dAs[0],cs[0],dcs[0],sprime[0],\
dsprime[0]])
print('FOM %=',FOM)
print(myTable)