# Deconvolution with transformed GOK-LM equation plus bgd
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
from prettytable import PrettyTable
import warnings
warnings.filterwarnings("ignore")
data = np.loadtxt('K120.txt')
x_data,y_data = data[:, 0], data[:, 1] 
y_data=y_data/max(y_data)
def GOKLM(t, Imax,b,tmax): 
    LM=Imax*t/tmax*(((b-1)/b)*(t**2/(tmax**2))/2+(b+1)/\
(2*b))**(b/(1-b))
    LM=LM+bgd*t/P
    return LM
P=int(max(x_data))
t = np.linspace(0, P, P)
bgd=y_data[-1]
inis=[Imax,b,tmax]=[1, 1.5, 60]
params, cov = optimize.curve_fit(GOKLM,\
x_data,y_data,p0=inis,maxfev=10000)   
res=GOKLM(x_data, *params)-y_data
FOM=round(100*np.sum(abs(res))/np.sum(y_data),2)
print('FOM=',FOM,' %')
[Imaxs,cs,tmaxs]= [round(params[x],3) for x in range(3)]
[dImaxs,dcs,dtmaxs]= [round(np.sqrt(cov[x][x]),3) for x in \
range(3)]
myTable = PrettyTable([ "Im (a.u.)","dIm",\
'b','db','tm (s)','dtm']) 
myTable.add_row([Imaxs,dImaxs,cs,dcs,tmaxs,dtmaxs])
print(myTable)