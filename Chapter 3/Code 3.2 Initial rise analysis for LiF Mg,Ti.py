# Initial rise analysis for LiF:Mg,Ti
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
from prettytable import PrettyTable 
import warnings
warnings.filterwarnings("ignore")
from scipy import stats
import os
os.chdir('C:/Users/Bill Pagonis/Desktop/pythonVP')
files=('LiFTL150degC.txt','LiFcurve3.txt','LiFcurve5.txt')
plt.subplot(1,2, 1)
plt.title('(a)')
plt.ylabel('TL signal [a.u.]')
plt.xlabel(r'Temperature [$^{o}$]')
sym=(['o-','s-','^-']) 
sym2=(['o','s','^'])
lab=('150$^{o}$C','164$^{o}$C','178$^{o}$C')
for i in range(3):
    data = np.loadtxt(files[i])
    x_data,y_data = data[:, 0], data[:, 1]
    plt.plot(x_data,y_data,sym[i],label=lab[i],linewidth=3)
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.subplot(1,2, 2)    
plt.title('(b)  Initial rise')
plt.ylabel('TL [a.u.]')
plt.xlabel(r'Temperature T [$^{o}$C]')
Evalues=[0 for x in range(len(files))]
Rvalues=[0 for x in range(len(files))]
dEvalues=[0 for x in range(len(files))]
lowpts=([3,3,3])
hipts=([9,9,9])
for i in range(3):
    data = np.loadtxt(files[i])
    x_data,y_data = data[:, 0][lowpts[i]:hipts[i]],\
			data[:, 1][lowpts[i]:hipts[i]]
    x_data=[1/(8.617e-5*(273+u)) for u in np.array(x_data)]
    y_data=np.log(y_data)
    plt.plot(x_data,y_data,sym2[i])
    slope, intercept,r_value,p_value, std_err=\
			stats.linregress(x_data,y_data)
    fit=[x*slope+intercept for x in x_data]
    plt.plot(x_data,fit,'-',label=lab[i])
    slope,std_err= np.round(slope,2), np.round(std_err,2)  
    Evalues[i]=round(slope,3)
    dEvalues[i]=round(std_err,3)
avgE=round(np.mean(Evalues),3)
dE=round(np.std(Evalues),2)
plt.tight_layout()
myTable = PrettyTable([ "Curve","E (eV)","dE (eV)",\
"average E (eV)","stdE (eV)"]) 
myTable.add_row(["1",Evalues[0],dEvalues[0],avgE,dE]);
myTable.add_row(["2",Evalues[1],dEvalues[1],' ',' ']);
myTable.add_row(["3",Evalues[2],dEvalues[2],' ',' ']);
print(myTable)
plt.show()