# superlinear dose response  of anion deficient aluminum oxide
from scipy.special import lambertw
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from prettytable import PrettyTable 
## fit to Lambert equation  ----
## fit to saturation exponential ----
t = ([0.0537568,0.103385,0.156481,0.211929,0.260776, 
0.321015,0.36483,0.414625,0.478926,0.535569,0.589476,
0.638899,0.824344,0.951985,1.18991,1.63688,3.19332])
y = ([6694.89,15592.6,24767.6,39360.5,52176.2,78075.6,
101463,131855,189548,227223,272406,376197,469268,
634929,792121,1.16105e6,1.6992e6])
x_data=np.array(t)
y_data=np.array(y)
#plt.plot(x_data,y_data)
def lambertfit(x_data,N,B,Dc,beta):
    u=np.real(lambertw((np.abs(B))*np.exp(np.abs(B)-(x_data/\
    np.abs(Dc))))/(np.abs(B)))
    u=N*(1-(u**beta))
    u.astype(float)
    return u
init_vals=[2e6,10,.01,.01]
params, cov = optimize.curve_fit(lambertfit,\
x_data, y_data,p0=init_vals)
dN = int(np.sqrt(cov[0][0]))
dB = round(np.sqrt(cov[1][1]),2)
dDc = round(np.sqrt(cov[2][2]),3)
dbeta = round(np.sqrt(cov[2][2]),2)
x_vals=np.arange(0.03,4,.01)
plt.subplot(1,2, 1)
plt.plot(x_vals, lambertfit(x_vals, *params[0:5]),c="b",\
label='Lambert fit') 
plt.title('(a)')   
plt.scatter(x_data, y_data, label='Experiment')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('OSL (L/T)')
plt.xlabel('Dose [Gy]')
plt.text(1.8,.7e6,'Anion-deficient')
plt.text(1.8,.6e6,'Aluminum Oxide')
plt.text(1.8,.5e6,'Linear scale')
plt.subplot(1,2, 2)
plt.title('(b)')
plt.plot(np.log(x_vals), np.log(lambertfit(x_vals, \
*params[0:5])),c="b",label='Lambert fit')     
plt.scatter(np.log(x_data), np.log(y_data), label='Experiment')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylabel('OSL (L/T)')
plt.xlabel('Dose [Gy]')
plt.text(-1,10,'Log-log scale')
plt.title('(b)')
plt.tight_layout()
myTable = PrettyTable(["N", "B","dB", "Dc (Gy)",\
"d(Dc","beta","dbeta"]) 
myTable.add_row([format(params[0],"10.1E"),\
round(np.abs(params[1]),2),dB,round(np.abs(params[2])\
, 2),dDc,round(np.abs(params[3]),5),dbeta])
print(myTable)
plt.show()