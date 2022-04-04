# Deconvolution of single TL peak with Lambert-OTOR equation
from scipy import optimize
import numpy as np
from prettytable import PrettyTable 
import warnings
warnings.filterwarnings("ignore")
from scipy.special import wrightomega
data = np.loadtxt('glocanin1.txt')
x_data,y_data = data[:, 0], data[:, 1]/max( data[:, 1])
kB, beta=  8.617e-5,  1
def W_func(T,A,sprime, E, c): 
   expint=kB*(T**2.0)/(beta*E)*np.exp(-E/(kB*T))*(1-2*kB*T/E)
   zTL=(1/c)-np.log(c)+(sprime*expint)
   lam=wrightomega(zTL)
   return A*np.exp(-E/(kB*T))/(lam+lam**2)                               
params, cov=optimize.curve_fit(W_func,x_data,y_data,([1e10,2e9,\
1.,10]),bounds=((1e8,1e8,.9,1e-4),(1e15,1e14,1.3,1e4)))
A=format(params[0],"10.1E")
sprime=format(params[1],"10.1E")
E=round(params[2],3)
c=round(params[3],1)
dA = format(np.sqrt(cov[0][0]),"10.1E")
dsprime = format(np.sqrt(cov[1][1]),"10.1E")
dE = round(np.sqrt(cov[2][2]),5)
dc = round(np.sqrt(cov[3][3]),1)
res=W_func(x_data, *params)-y_data
FOM=round(100*np.sum(abs(res))/np.sum(y_data),1)
myTable = PrettyTable()
myTable.add_row([ "A","dA", "s' (s^-1)","ds' (s^-1)"])
myTable.add_row([A,dA,sprime,dsprime])
myTable.add_row([" "]*4)
myTable.add_row(["E (eV)","dE(eV)","c","dc"])
myTable.add_row([E,dE,c,dc])
print("FOM=",FOM," %")
print(myTable)