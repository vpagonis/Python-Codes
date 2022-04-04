#The effect of thermal quenching on TL glow curve for quartz
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
from scipy import stats
def fn1():
    for i in range(len(files)):
        data = np.loadtxt(files[i])
        x_data,y_data = data[:, 0], data[:, 1]
        plt.plot(x_data,y_data,'o-')
        area[i]=np.sum([x_data[u]*y_data[u] for u in \
		range(len(x_data))])
        Tmax[i]=x_data[np.argmax(y_data)]
def fn2():
    for i in range(len(files)):
        data = np.loadtxt(files[i])
        x_data,y_data = data[:, 0][lowpts[i]:hipts[i]],\
			data[:, 1][lowpts[i]:hipts[i]]
        x_data=[1/(8.617e-5*u) for u in np.array(x_data)]
        y_data=np.log(y_data)
        plt.plot(x_data,y_data)
        slope, intercept,r_value,p_value, std_err=\
			stats.linregress(x_data,y_data)
        fit=[x*slope+intercept for x in x_data]
        plt.plot(x_data,fit,'o-')
        slope,std_err= np.round(slope,2), np.round(std_err,2)  
        Evalues[i]=round(slope,3)
        dEvalues[i]=round(std_err,3)
files=('B2qz025.txt','B2qz05.txt','B2qz2.txt','B2qz5.txt',\
'B2qz10.txt')
plt.subplot(221)
plt.title('(a)')
plt.ylabel('TL signal [a.u.]')
plt.xlabel(r'Temperature [K]')
plt.ylim(0,22)
plt.text(350,17,'0.25-10 K/s')
area=[0 for x in range(len(files))]
Tmax=[0 for x in range(len(files))]
fn1()
plt.subplot(222)    
plt.title('(b)  Initial rise')
plt.ylabel('TL [a.u.]')
plt.xlabel(r'Temperature [K]')
Evalues=[0 for x in range(len(files))]
Rvalues=[0 for x in range(len(files))]
dEvalues=[0 for x in range(len(files))]
lowpts=[3 for x in range(len(files))]
hipts=[9 for x in range(len(files))]
fn2()
plt.subplot(223)
plt.plot(Tmax,np.array(area)/1e5,'o')
plt.xlabel(r'T$_{max}$ [K]')
plt.ylabel(r'Area of TL peak [a.u.]')
plt.text(330,.2,"Thermal quenching")
plt.text(330,.08,"of TL area")
plt.title('(c)')
plt.ylim(0,1)
plt.subplot(224)
plt.ylim(0,1.2)
plt.title('(d)')
plt.text(330,.3,"Thermal quenching")
plt.text(330,.15,"of Energy E")
plt.plot(Tmax,-np.array(Evalues),'r^')
plt.xlabel(r'T$_{max}$ [K]')
plt.ylabel(r'Energy E [eV]')
plt.tight_layout()
plt.show()