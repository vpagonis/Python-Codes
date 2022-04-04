# Deconvolution of Glocanin TL with original GOK-TL
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
from scipy import optimize
from prettytable import PrettyTable 
data = np.loadtxt('glocanin1.txt')
x_data,y_data = data[:, 0], data[:, 1]/max(data[:, 1])
kB,  beta= 8.617e-5,  1
def TLGOK(T,A,sprime,E,b):    
    expint=kB*(T**2.0)/(beta*E)*\
	    np.exp(-E/(kB*T))*(1-2*kB*T)/E
    return A*np.exp(-E/(kB*T))*((1+sprime*(b-1)*expint)\
    **(-b/(b-1)) )
params, cov=optimize.curve_fit(TLGOK,x_data,y_data,([1,2e10,1.,\
1.5]),bounds=((1e-10,1e8,1,1),(1e20,1e14,1.3,1.999)))
A=format(params[0],"10.1E")
dA = format(np.sqrt(cov[0][0]),"10.1E")
sprime=format(params[1],"10.1E")
E=round(params[2],4)
dE = round(np.sqrt(cov[2][2]),7)
b=round(params[3],5)
db = round(np.sqrt(cov[3][3]),7)
res=TLGOK(x_data, *params)-y_data
FOM=round(100*np.sum(abs(res))/np.sum(y_data),2)
myTable = PrettyTable([ "A","s' (s^-1)","E(eV)",\
"dE","b","db"]) 
myTable.add_row([A,sprime,E,dE,b,db]);
print("FOM=",FOM," %")
print(myTable)