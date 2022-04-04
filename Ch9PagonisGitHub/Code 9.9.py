# Deconvolution with MOK-LM equation, fixed bgd
from scipy import optimize
import numpy as np
from prettytable import PrettyTable
import warnings
warnings.filterwarnings("ignore")
data = np.loadtxt('K120.txt')
x_data,y_data = data[:, 0], data[:, 1] 
y_data=y_data/max(y_data)
def MOKLM(t, A,b,sprime): 
    F=np.exp(sprime*t**2/(2*P))
    LM=A*t/P*F/((F-b)**2) 
    LM=LM+bgd*t/P
    return LM
P=int(max(x_data))
t = np.linspace(0, P, P)
inis=[A,b,sprime]=[.5, .1, 1e-4]
bgd=y_data[-1]
params, cov = optimize.curve_fit(MOKLM,\
x_data,y_data,p0=inis,maxfev=10000)  
res=MOKLM(x_data, *params)-y_data
FOM=round(100*np.sum(abs(res))/np.sum(y_data),2)
res=MOKLM(x_data, *params)-y_data
FOM=round(100*np.sum(abs(res))/np.sum(y_data),2)
print('FOM=',FOM,' %')
[As,cs,sprimes]= [round(params[x],3) for x in range(3)]
[dAs,dcs,dsprimes]= [round(np.sqrt(cov[x][x]),3)\
for x in range(3)]
myTable = PrettyTable([ "A (a.u.)","dA",\
'alpha','dalpha','s (s^-1)','ds']) 
myTable.add_row([As,dAs,cs,dcs,sprimes,dsprimes])
print(myTable)