#Fit dose response data with Saturating Exponential
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import lambertw
import warnings
warnings.filterwarnings("ignore")
from scipy import optimize
from prettytable import PrettyTable 
## fit to SE equation  ----
x_data= ([0, 3.55, 44.18, 258.718, 1051.62, 2044.98, 3003.94,
      5024.61, 7046.32, 9992.29])
y_data = ([0, 0.93512, 2.61108, 4.99104, 6.36704, 6.42148, 
      6.43643, 6.46792, 6.77215, 6.97391])
# x_data =np.array([0,50.7117,100.534,152.135,204.626,272.242])
# y_data = np.array([0,33.144,42.205,43.1055,44.4157,43.7098])
plt.plot(x_data,y_data,"o")
def SE(x_data,A,Do):
    u=abs(A)*(1-np.exp(-x_data/np.abs(Do)))
    return u
def lambertfit(x_data,N,R,Dc):
    u=np.real(N*(1+lambertw((np.abs(R)-1)*np.exp(np.abs(R)-1-\
    x_data/np.abs(Dc)))/(1-np.abs(R))))
    u.astype(float)
    return u
def DSE(x_data,A,Do1,B,Do2):
    u=A*(1-np.exp(-x_data/np.abs(Do1)))+B*(1-\
    np.exp(-x_data/np.abs(Do2)))
    return u
init_vals=[1,1000,1,10000]                #DSE
params, cov = optimize.curve_fit(DSE,\
x_data, y_data,p0=init_vals)
[A,Do]=[round(params[x],1) for x in range(2)]
dA = round(np.sqrt(cov[0][0]),2)
dDo = int(np.sqrt(cov[1][1]))
x_vals=np.arange(0,1e4,1e1)
plt.plot(x_vals, DSE(x_vals, *params[0:4]),c="b",\
linestyle="dashed",label='DSE fit') 
plt.xlim(0,3200)
[A1,Do1,A2,Do2]=[abs(round(params[x],1)) for x in range(4)]
dA1 = round(np.sqrt(cov[0][0]),2)
dDo1 = int(np.sqrt(cov[1][1]))
dA2 = round(np.sqrt(cov[2][2]),2)
dDo2 = int(np.sqrt(cov[3][3]))
myTable = PrettyTable([" ","A1", "dA1", "Do1","dDo1","A2",\
"dA2","D2","dD2"]) 
myTable.add_row(["DSE",A1,dA1,Do1,\
dDo1,A2,dA2,Do2,dDo2])
print(myTable)
init_vals=[20,20]                        #SE
params, cov = optimize.curve_fit(SE,\
x_data, y_data,p0=init_vals)
[A,Do]=[round(params[x],1) for x in range(2)]
dA = round(np.sqrt(cov[0][0]),2)
dDo = int(np.sqrt(cov[1][1]))
plt.plot(x_vals, SE(x_vals, *params[0:2]),c="r",\
linestyle="dotted",label='SE fit')     
plt.scatter(x_data, y_data, label='Experiment')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('OSL (L/T)')
plt.xlabel('Dose [Gy]')
plt.title('Coarse grain quartz')
dA = round(np.sqrt(cov[0][0]),1)
dDo = round(np.sqrt(cov[1][1]),1)
init_vals=[100.0,.1,100]
params, cov = optimize.curve_fit(lambertfit,\
x_data, y_data,p0=init_vals)
dN = round(np.sqrt(cov[0][0]),0)
dR = format(abs(np.sqrt(cov[1][1])),"10.1E")
dDc = int(np.sqrt(cov[2][2]))
plt.plot(x_vals, lambertfit(x_vals, *params[0:4]),\
label='PKC fit')     
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.tight_layout()
myTable2 = PrettyTable([" ","N", "dN", "Do (Gy)","dDc"]) 
myTable2.add_row(["SE",A,dA,Do,\
dDo])
print(myTable2)
myTable3=PrettyTable([" ","N","dN", "R","dR","Dc (Gy)","dDc"]) 
myTable3.add_row(["PKC",round(params[0],0),dN,\
format(abs(params[1]),"10.1E"),dR,abs(int(params[2])),dDc])
print(myTable3)
plt.show()