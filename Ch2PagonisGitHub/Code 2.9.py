# Apply the initial rise method to find the activation energy E
# Load the data from txt file, which  contains pairs of
# data  in the form: (Temperature_in_K,TL_Intensity (any units)
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
from prettytable import PrettyTable 
data = np.loadtxt('lbodata.txt')
x_data,y_data = data[:, 0], data[:, 1]
plt.subplot(1,2, 1)
plt.plot(x_data,y_data,'+-',c='green',label='TL')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('TL signal [a.u.]')
plt.xlabel(r'Temperature T [K]')
plt.ylim(0,140)
plt.text(370, 100,'Sample')
plt.text(370, 90,'LBO')
kB=8.617*1e-5 # Boltzmann constant in eV/K
initialPos=7  #analyze data points  #6 to #14
finalPos=13
plt.plot(x_data[initialPos:finalPos],y_data[initialPos:\
finalPos],"o",c='r',label='IR-data')
plt.title('(a)')
x= 1/(kB*x_data[initialPos:finalPos])
y= np.log(y_data[initialPos:finalPos])
slope, intercept,r_value,p_value, std_err=stats.linregress(x,y)
plt.subplot(1,2,2)
plt.plot(x,y,'ro',label='Data')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.plot(x,x*slope+intercept,label='Fit')
plt.xlabel(r'1/(kT) (eV$^{-1}$)')
plt.ylabel('ln(TL)')
plt.title('(b)  Initial rise method')
plt.tight_layout()
slope, intercept, rsqr, p_value,std_err=\
np.round(slope,2),np.round(intercept,2),np.round(\
r_value**2.0,3),format(p_value,"10.2E"), np.round(std_err,2)
myTable=PrettyTable(['E (eV)','dE (eV)','p-value','R^2']) 
myTable.add_row([-slope, std_err, p_value,rsqr])
print(myTable)
plt.show()