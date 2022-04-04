#Fit dose response data with Saturating Exponential
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import lambertw
import warnings
warnings.filterwarnings("ignore")
from scipy import optimize
from prettytable import PrettyTable 
## fit to SE equation  ----
x_data = ([-34.2, 34.2466,68.4932,273.973, 1027.4,1986.3,3013.7,
      5000, 7979.45, 10000])
y_data = ([1.04, 0.000978474, 2.76386, 7.24592,12.6008,14.5329,
      15.8956,  17.1905, 17.847, 18.0952])
# x_data =np.array([0,50.7117,100.534,152.135,204.626,272.242])
# y_data = np.array([0,33.144,42.205,43.1055,44.4157,43.7098])
plt.plot(x_data,y_data,"o")
def SE(x_data,A,Do):
    u=A*(1-np.exp(-x_data/np.abs(Do)))
    return u
def lambertfit(x_data,N,R,Dc):
    u=np.real(N*(1+lambertw((np.abs(R)-1)*np.exp(np.abs(R)-1-\
    x_data/np.abs(Dc)))/(1-np.abs(R))))
    u.astype(float)
    return u
init_vals=[20,20]
params, cov = optimize.curve_fit(SE,\
x_data, y_data,p0=init_vals)
[A,Do]=[round(params[x],1) for x in range(2)]
dA = round(np.sqrt(cov[0][0]),2)
dDo = int(np.sqrt(cov[1][1]))
x_vals=np.arange(0,1e4,1e2)
plt.plot(x_vals, SE(x_vals, *params[0:2]),linestyle='--',c="b",\
label='SE fit')     
plt.scatter(x_data, y_data, label='Experiment')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylim(0,30)
plt.ylabel('OSL (L/T)')
plt.xlabel('Dose [Gy]')
plt.title('Fine grain dose response')
dA = round(np.sqrt(cov[0][0]),1)
dDo = round(np.sqrt(cov[1][1]),1)
init_vals=[100.0,.1,100]
params, cov = optimize.curve_fit(lambertfit,\
x_data, y_data,p0=init_vals)
dN = round(np.sqrt(cov[0][0]),0)
dR = round(abs(np.sqrt(cov[1][1])),1)
dDc = int(np.sqrt(cov[2][2]))
plt.plot(x_vals, lambertfit(x_vals, *params[0:4]),\
label='PKC fit')     
plt.scatter(x_data, y_data)
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.tight_layout()
myTable = PrettyTable(["-","N", "dN", "R","dR","Do/Dc (Gy)","dDc"]) 
myTable.add_row(["SE",A,dA,' ','',Do,\
dDo])
myTable.add_row(["PKC",round(params[0],0),dN,\
round(abs(params[1]),1),dR,abs(int(params[2])),dDc])
print(myTable)
plt.show()