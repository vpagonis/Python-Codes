## Fit of quartz OSL dose response data using W(x)
from scipy.special import lambertw
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from prettytable import PrettyTable 
t = ([-34.2466, 34.2466, 68.4932, 273.973, 1027.4,1986.3,3013.7,
      5000, 7979.45, 10000])
y = ([1.04664, 0.000978474, 2.76386, 7.24592,12.6008,14.5329,
      15.8956,  17.1905, 17.847, 18.0952])
t2= ([0, 3.5583, 44.1822, 258.718, 1051.62, 2044.98, 3003.94,
      5024.61, 7046.32, 9992.29])
y2 = ([0, 0.93512, 2.61108, 4.99104, 6.36704, 6.42148, 
      6.43643, 6.46792, 6.77215, 6.97391])
x_data=np.array(t)
y_data=np.array(y)
#plt.plot(x_data,y_data)
def lambertfit(x_data,N,R,Dc):
    u=np.real(N*(1+lambertw((np.abs(R)-1)*np.exp(np.abs(R)-1-\
    x_data/np.abs(Dc)))/(1-np.abs(R))))
    u.astype(float)
    return u
init_vals=[20.0,.1,1e4]
params, cov = optimize.curve_fit(lambertfit,\
x_data, y_data,p0=init_vals)
dN = round(np.sqrt(cov[0][0]),2)
dR = format(np.sqrt(cov[1][1]),"10.2E")
dDo = round(np.sqrt(cov[2][2]),2)
x_vals=np.arange(0,1e4,1)
plt.subplot(1,2, 1)
plt.plot(x_vals, lambertfit(x_vals, *params[0:4]),\
label='Lambert fit')     
plt.scatter(x_data, y_data, label='Experiment')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('OSL (L/T)')
plt.xlabel('Dose [Gy]')
plt.text(3000,10,'Fine grain quartz')
plt.title('(a)')
plt.tight_layout()
myTable=PrettyTable(["N","dN","R","dR","Dc(Gy)","d(Dc)(Gy)"]) 
myTable.add_row([round(params[0],1),dN,format(np.abs(\
params[1]),"10.1E"),dR,np.int(np.abs(params[2])),np.int(dDo)])
x_data=np.array(t2)
y_data=np.array(y2)
init_vals=[20.0,.1,1e4]
params, cov = optimize.curve_fit(lambertfit,\
x_data, y_data,p0=init_vals)
dN = round(np.sqrt(cov[0][0]),2)
dR = format(np.sqrt(cov[1][1]),"10.2E")
dDo = round(np.sqrt(cov[2][2]),2)
x_vals=np.arange(0,1e4,1)
plt.subplot(1,2, 2)
plt.plot(x_vals, lambertfit(x_vals, *params[0:4]),\
label='Lambert fit')     
plt.scatter(x_data, y_data, label='Experiment')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('OSL (L/T)')
plt.xlabel('Dose [Gy]')
plt.text(2500,4,'Coarse grain quartz')
plt.title('(b)')
plt.tight_layout()
myTable.add_row([round(params[0],1),dN,format(np.abs(\
params[1]),"10.1E"),dR,np.int(np.abs(params[2])),np.int(dDo)])
print(myTable)
plt.show()