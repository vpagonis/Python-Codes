# J1000 deconvolution with optimal number of peaks N=7 
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
from prettytable import PrettyTable 
import warnings
warnings.filterwarnings("ignore")
data = np.loadtxt('J1000prompt.txt')
x_data,y_data = data[:, 0], data[:, 1] 
z,   kB =1.8, 8.617E-5
def TL(T, B,En ,rho, s):
    return abs(B)* np.exp(-rho*( (np.log(1+z*s*kB*(((T+\
				    273)**2.0)/np.abs(En))*np.exp(-En/(kB*(T+273)))*\
				    (1-2*kB*(T+273)/En)))**3.0))*(En**2.0-6*(kB**2.0)*\
    	((T+273)**2.0))*((np.log(1+z*s*kB*(((T+273)**2.0)/\
				    abs(En))*np.exp(-En/(kB*(T+273)))*(1-2*kB*(T+273)/\
				    En)))**2.0)/(En*kB*s*((T+273)**2)*z-2*(kB**2.0)*\
				    s*z*((T+273)**3.0)+np.exp(En/(kB*(T+273)))*En)
def total_TL(T, *inis): 
    u=np.array([0 for i in range(len(x_data))])
    Bs, Ens,rho,s=    inis[0:nPks], inis[nPks:2*nPks],inis[-2]\
		    ,inis[-1]
    for i in range(nPks):
        u=u+TL(T,Bs[i],Ens[i],rho,s)
    return u
nPks=7
B=[2e17]*7
lowB, highB=[0.01*x for x in B], [50*x for x in B]
lowrho, highrho, rho= [0.003,.015,.008]
lows, highs, s =[1e11, 1e14, 1e12]
En=[.8, .9,1.06,1.19,1.4,1.6,1.7]
lowEn,highEn =[0.9*x for x in En],[1.1*x for x in En]
inis=B+En+[rho]+[s]
lowbnds=lowB+lowEn+[lowrho]+[lows]
highbnds=highB+highEn+[highrho]+[highs] 
params, params_covariance = optimize.curve_fit(total_TL,\
x_data,y_data,p0=inis,bounds=(lowbnds,highbnds),maxfev=10000)
plt.scatter(x_data, y_data,c='r',label='Sample J1000')
plt.plot(x_data, total_TL(x_data, 
 *params),c='black',label='KP-TL  N='+str(nPks),linewidth=1)   
for i in range(0,nPks): 
    plt.plot(x_data, TL(x_data, params[i],params[nPks+i],\
		    params[-2],params[-1]))
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)  
plt.ylabel('TL [a.u.]')
plt.xlabel(r'Temperature T [$^{o}$C]')
res=total_TL(x_data, *params)-y_data
FOM=100*np.sum(abs(res))/np.sum(y_data)
plt.text(90,.2e6,'1')
plt.text(110,1e6,'2')
plt.text(150,1e6,'3')
plt.text(210,.4e6,'4')
plt.text(280,.3e6,'5')
plt.text(330,.2e6,'6')
plt.text(390,.2e6,'7')
plt.show()