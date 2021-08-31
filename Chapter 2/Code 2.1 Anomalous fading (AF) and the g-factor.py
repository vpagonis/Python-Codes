# Anomalous fading (AF) and calculating  the g-factor
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from prettytable import PrettyTable 
def linearFunc(x2,intercept,slope):
    y2 = intercept + slope * x2
    return y2
import os
os.chdir('C:/Users/Bill Pagonis/Desktop/pythonVP')
import os
os.chdir('C:/Users/Bill Pagonis/Desktop/pythonVP')
data = np.loadtxt('durnago0sok.txt')
x0, y0 = np.array(data[:, 0]), np.array(data[:, 1])
data = np.loadtxt('durango10daysok.txt')
x1, y1 = np.array(data[:, 0]), np.array(data[:, 1])
data = np.loadtxt('DurangoAFdataok.txt')
x2, y2 = np.array(data[:, 0]), np.array(data[:, 1])
plt.subplot(1,2, 1)
plt.plot(x0,y0,'b^-',label='t=0s')
plt.plot(x1,y1,'ro-',label='t=10 days')
plt.title('(a)')
plt.ylim(0,1.1e5);
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('Remnant TL signal [a.u.]')
plt.xlabel(r'Temperature [$^{o}$C]')
plt.text(100, 6e4,'Durango')
plt.text(100,5e4,'apatite')
#######
plt.subplot(1,2, 2)
plt.plot(x2,y2,"o",c="r")
plt.ylabel('Remnant TL signal [a.u.]')
plt.xlabel(r'ln(time)')
plt.title('(b)')
plt.text(5,.9,"g-factor")
plt.text(5,.85,r"analysis")
plt.text(5,.76,r"g=21.0%")
params,cov=curve_fit(linearFunc,x2,y2)
plt.plot(x2,linearFunc(x2,*params),label='Fit')
inter,slope = np.round(params,4)
d_inter,d_slope = np.round(np.sqrt(np.diag(cov)),3)
residuals = y2- linearFunc(x2, *params)
ss1=np.sum(residuals**2.0)
ss2 = np.sum((y2-np.mean(y2))**2.0)
rsqr=np.round(1-ss1/ss2,3)
g=np.round(-230.2*slope)
myTable=PrettyTable(['g (%)','slope (a.u.)','dslope (a.u.)',\
'R^2']) 
myTable.add_row([g,-slope, d_slope, rsqr]);
print(myTable)
plt.tight_layout()
plt.show()