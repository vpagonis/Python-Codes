# Fit to PKC-S equation for Supralinearity index f(D)
from scipy.special import lambertw
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
from scipy import optimize
from prettytable import PrettyTable 
## fit to PKC equation  ----
t = ([0.0811131, 0.171804,0.450923,0.857988,1.88341,      
4.44069,8.44947,17.2683,40.7152,83.2104,176.246,
415.552,819.456,1674.74,3948.68,8983.32])
y = ([1.03326,0.983911,1.08074,1.10465,1.12844,1.28023,
1.37731,1.50481,1.76636,1.79021,1.45428, 0.813382,      
0.428722,0.208669,0.0799713,0.0305689])
x_data=np.array(t)
y_data=np.array(y)
def lambertfit(x_data,N,B,Dc,beta):
    u=np.real(lambertw((np.abs(B))*np.exp(np.abs(B)-(x_data/\
np.abs(Dc))))/(np.abs(B)))
    u=N*(1-(u**beta))/x_data
    u.astype(float)
    return u
init_vals=[100,10,1,.01]
params, params_covariance = optimize.curve_fit(lambertfit,\
x_data, y_data,p0=init_vals)
x_vals=np.arange(0.1,1e4,.1)
plt.subplot(1,2, 1)
plt.plot(x_vals,lambertfit(x_vals,*params[0:5]),c="b",\
label='Lambert fit')     
plt.scatter(x_data, y_data, label='Experiment')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('Supralinearity index f(D)')
plt.title('(a)')
plt.xlabel('Dose [Gy]')
plt.text(4000,.65,'Probe A')
plt.text(4000,.5,'Linear scale')
plt.subplot(1,2, 2)
plt.title('(b)')
plt.plot(np.log(x_vals),lambertfit(x_vals,*params[0:5]),\
c="b",label='Lambert fit')     
plt.scatter(np.log(x_data), y_data, label='Experiment')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('Supralinearity index f(D)')
plt.xlabel('ln(Dose) [KGy]')
plt.text(-2,.5,'Log-Linear scale')
plt.title('(a)')
plt.tight_layout()
res=lambertfit(x_data,*params[0:5])-y_data
FOM=round(100*np.sum(abs(res))/np.sum(y_data),2)
myTable = PrettyTable(["N", "B", "Dc (KGy)","beta","FOM (%)"]) 
myTable.add_row([round(params[0],2),format(np.abs(params[1]),\
"10.2E"),format(np.abs(params[2]), "10.2E"),round(np.abs(\
params[3]),4),FOM])
print(myTable)
plt.show()