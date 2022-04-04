#  Analysis of CW-IRSL measured at high stimulation temeprature 
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
[betas, dbetas, As,dAs,taus,dtaus,maxys,bgds,dbgds]= [[] \
for _ in range(9)]
fils=['irrthenst50.txt','irrthenst75.txt','irrthenst100.txt',
'irrthenst125.txt','irrthenst150.txt','irrthenst175.txt',
'irrthenst200.txt','irrthenst225.txt'] 
def total_CW(t, *inis): 
    [As, taus, betas,bgd]=  inis
    CW=As*np.exp(-(t/abs(taus))**np.abs(betas))+np.abs(bgd)
    return CW
def findbeta(fil):
    data = np.loadtxt(fil)
    y_data=data
    x_data=np.arange(0,len(y_data))
    maxys.append(max(y_data))
    y_data=y_data/max(y_data)
    inis=[1,5,.1,min(y_data)]
    params, cov = optimize.curve_fit(total_CW,\
    x_data,y_data,p0=inis,maxfev=10000)   
    res=total_CW(x_data, *params)-y_data
    FOM=100*np.sum(abs(res))/np.sum(y_data)
    As.append(params[0])
    taus.append(params[1])
    betas.append(params[2])
    bgds.append(params[3])
    dAs.append(np.sqrt(cov[0][0]))
    dtaus.append(np.sqrt(cov[1][1]))
    dbetas.append(np.sqrt(cov[2][2]))
    dbgds.append(np.sqrt(cov[3][3]))
    return
for j in range(len(fils)): 
    findbeta(fils[j])
temps=[50,75,100,125,150,175,200,225]
plt.subplot(2,2,1)
plt.errorbar(temps, maxys, marker='s',yerr=np.sqrt(maxys),\
label=r'Amplitude A$')
plt.ylim(0,6e5)
plt.ylabel(r'Amplitude A$')
plt.xlabel(r'Stimulation T [$^o$C]')
plt.subplot(2,2,2)
plt.errorbar(temps, betas, marker='s',yerr=dbetas,\
label=r'Coefficient $\beta$')
plt.ylim(0.6,.7)
plt.ylabel(r'Coefficient $\beta$')
plt.xlabel(r'Stimulation T [$^o$C]')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0) 
plt.subplot(2,2,3)
plt.errorbar(temps, np.abs(taus), marker='s',yerr=dtaus)
plt.ylim(0,8)
plt.ylabel(r'$\tau$ [s]')
plt.xlabel(r'Stimulation T [$^o$C]')
plt.subplot(2,2,4)
plt.errorbar(temps,np.abs(bgds),marker='o',yerr=dbgds)
plt.ylabel(r'bgd [a.u.]')
plt.xlabel(r'Stimulation T [$^o$C]')
plt.ylim(0,.022)
plt.tight_layout()
plt.show()