# Deconvolution with transformed MOK-LM equation, fixed bgd
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
from prettytable import PrettyTable
import warnings
warnings.filterwarnings("ignore")
data = np.loadtxt('K120.txt')
x_data,y_data = data[:, 0], data[:, 1] 
y_data=y_data/max(y_data)
# b in this code represents the MOK parameter alpha
def MOKLM(t, Imax,b,tmax): 
    Fm=1.6476-1.0012*b+0.357*b**2
    F=np.exp((t/tmax)**2*(Fm-b)/(2*(Fm+b)))
    LM=Imax*t/tmax*(((Fm-b)/(F-b))**2)*F/Fm
    LM=LM+bgd*t/P
    return LM
P=int(max(x_data))
t = np.linspace(0, P, P)
bgd=y_data[-1]
inis=[Imax,b,tmax]=[1, .8, 60]
params, cov = optimize.curve_fit(MOKLM,\
x_data,y_data,p0=inis,maxfev=10000)   
res=MOKLM(x_data, *params)-y_data
FOM=round(100*np.sum(abs(res))/np.sum(y_data),2)
print('FOM=',FOM,' %')
[Imaxs,cs,tmaxs]= [round(params[x],3) for x in range(3)]
[dImaxs,dcs,dtmaxs]= [round(np.sqrt(cov[x][x]),3) for x in \
range(3)]
myTable = PrettyTable([ "Im (a.u.)","dIm",\
'alpha','dalpha','tm (s)','dtm']) 
myTable.add_row([Imaxs,dImaxs,cs,dcs,tmaxs,dtmaxs])
print(myTable)