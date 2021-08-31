#Fit dose response data with Saturating Exponential
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from prettytable import PrettyTable 
## fit to SE equation  ----
x_data =np.array([0,50.7117,100.534,152.135,204.626,272.242])
y_data = np.array([0,33.144,42.205,43.1055,44.4157,43.7098])
plt.plot(x_data,y_data,"o")
def lambertfit(x_data,N,Do):
    u=N*(1-np.exp(-x_data/np.abs(Do)))
    return u
init_vals=[20,20]
params, cov = optimize.curve_fit(lambertfit,\
x_data, y_data,p0=init_vals)
dN = round(np.sqrt(cov[0][0]),2)
dDo = round(np.sqrt(cov[1][1]),2)
x_vals=np.arange(0,300,1)
plt.plot(x_vals, lambertfit(x_vals, *params[0:2]),c="b",\
label='SE fit')     
plt.scatter(x_data, y_data, label='Experiment')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylim(0,60);
plt.ylabel('OSL (L/T)')
plt.xlabel('Dose [Gy]')
plt.title('OSL dose response')
plt.text(200,20,'Libyan quartz')
plt.text(200,15,'Saturating exponential')
plt.tight_layout()
res=lambertfit(x_data, *params[0:2])-y_data
FOM=round(100*np.sum(abs(res))/np.sum(y_data),2)
myTable = PrettyTable(["N", "dN", "Do (Gy)","d(Do) (Gy)",\
"FOM (%)"]) 
myTable.add_row([round(params[0],2),dN,round(params[1],2),\
dDo,FOM])
print(myTable)
plt.show()