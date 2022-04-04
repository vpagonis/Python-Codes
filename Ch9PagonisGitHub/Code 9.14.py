#Dependence of CW-OSL signal on power 10-90% in Quartz
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
data = np.loadtxt('Quartz_110OSL_polymeris.txt')
lps=np.arange(1,10,1)
plt.subplot(1,2,1)
for i in lps:
    x_data, y_data = np.array(data[:, 0]), np.array(data[:, i])
    plt.plot(x_data,y_data)
plt.xlabel('Time [s]')
plt.ylabel('CW-IRSL intenisty I(t) [a.u.]')
plt.text(100,900,"CW-IRSL")
plt.text(100,800,'Variable LED power')
plt.subplot(1,2,2)    
for i in lps:
    x_data, y_data = np.array(data[:, 0]), np.array(data[:, i])
    y_data=np.array(x_data)*np.array(y_data)
    x_data =np.log(x_data)
    plt.plot(x_data,y_data)
plt.xlabel('ln(Time)')
plt.ylabel('t x I(t) [a.u.]')
plt.text(0,18000,"CW-IRSL")
plt.text(0,20000,'Trasformed')
plt.tight_layout()
plt.show()