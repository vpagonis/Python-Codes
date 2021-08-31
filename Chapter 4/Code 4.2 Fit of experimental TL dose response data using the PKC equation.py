from scipy.special import lambertw
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from prettytable import PrettyTable 
## fit to Lambert equation  ----
x_data =np.array([0.00394056,0.00451523,0.39225,0.41203,0.703,\
0.741318,0.742221,1.49553,1.49758,1.5158,2.98473,3.00304,\
3.02282, 5.99852, 6.05755])
y_data = np.array([2.45103,2.80847,3.97964,4.28586,5.30465,\
5.10007,5.66177,6.21706,7.49364,6.82965,8.50226,7.88934,\
8.19555,11.0809,11.7953])
x_data=np.array(x_data)
y_data=np.array(y_data)
plt.plot(x_data,y_data,"o")
def lambertfit(x_data,N,R,Dc,f):
    u=np.real(N*(1+lambertw((np.abs(R)-1)*np.exp(np.abs(R)-1-\
    (x_data+np.abs(f))/np.abs(Dc)))/(1-np.abs(R))))
    u.astype(float)
    return u
init_vals=[100.0,.1,1e6,2]
params, params_covariance = optimize.curve_fit(lambertfit,\
x_data, y_data,p0=init_vals)
x_vals=np.arange(-1,6,.1)
plt.plot(x_vals, lambertfit(x_vals, *params[0:5]),c="b",\
label='Lambert fit')     
plt.scatter(x_data, y_data, label='Experiment')
leg = plt.legend()
plt.ylim(0,13);
plt.xlim(-1,7);
leg.get_frame().set_linewidth(0.0)
plt.ylabel('OSL (L/T)')
plt.xlabel('Dose [Gy]')
plt.title('OSL dose response')
plt.text(3,4,'Volcanic glass')
plt.tight_layout()
res=lambertfit(x_data, *params[0:5])-y_data
FOM=round(100*np.sum(abs(res))/np.sum(y_data),2)
myTable = PrettyTable(["N", "R", "Dc (Gy)","f","FOM (%)"]) 
myTable.add_row([round(params[0],2),format(np.abs(params[1]),\
"10.2E"),format(np.abs(params[2]), "10.2E"),\
round(params[3],2),FOM])
print(myTable)
plt.show()