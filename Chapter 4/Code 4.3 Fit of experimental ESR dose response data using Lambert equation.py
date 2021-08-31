from scipy.special import lambertw
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from prettytable import PrettyTable 
## fit to Lambert equation  ----
TLqzx=([174.13, 345.027, 931.847, 1603.74, 2524.39, 4031.12,
6372.1,9044.1,12217.5,15058.3,19981.5,25072.3,30082.3,40011.4])
TLqzy=([1.07478, 1.68389, 2.18591, 2.93875, 3.51271, 4.48122,
5.59377,6.3484,7.3184,8.39557,9.33133,10.4105,11.8478,13.6478])
x_data=np.array(TLqzx)
y_data=np.array(TLqzy)
def lambertfit(x_data,N,R,Dc):
    u=np.real(N*(1+lambertw((np.abs(R)-1)*np.exp(np.abs(R)-1-\
    x_data/np.abs(Dc)))/(1-np.abs(R))))
    u.astype(float)
    return u
init_vals=[100.0,.1,1e6]
params, cov = optimize.curve_fit(lambertfit,\
x_data, y_data,p0=init_vals)
dN = round(np.sqrt(cov[0][0]),2)
dR = round(np.sqrt(cov[1][1]),2)
dDc = round(np.sqrt(cov[2][2]),2)
plt.plot(x_data, lambertfit(x_data, *params[0:4]),\
label='Lambert fit')     
plt.scatter(x_data, y_data, label='Experiment')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('ESR signal [a.u.]')
plt.xlabel('Dose [Gy]')
plt.title('Quartz ESR dose response')
plt.tight_layout()
res=lambertfit(x_data, *params[0:4])-y_data
FOM=round(100*np.sum(abs(res))/np.sum(y_data),1)
myTable=PrettyTable(["N","R","dR","Dc(Gy)",\
"d(Dc)(Gy)","FOM"]) 
myTable.add_row([round(params[0],1),format(np.abs(params[1]),\
"10.1E"),dR,format(np.abs(params[2]), "10.1E"),dDc,FOM])
print(myTable)
plt.show()