import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from prettytable import PrettyTable
fils=['OSL0s.txt','OSL100s.txt','OSL1000s.txt']
plt.subplot(1,2,1)
marks=['+','o','^']
labls=['0 s','100 s','1000 s']
for j in range(0,3):
    data = np.loadtxt(fils[j])
    x_data, y_data = 1+np.array(data[:, 0]), \
    np.array(data[:, 1])
    plt.plot(x_data,y_data,marks[j]+'-',label=labls[j])
plt.xlabel('Stimulation time t [s]')
plt.ylabel(r'I(t) CW-IRSL intensity  [cts/s]')
plt.title('(a)')
plt.text(500,2.5e6,'CW-IRSL signal')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0) 
plt.subplot(1,2,2)
for j in range(0,3):
    data = np.loadtxt(fils[j])
    x_data,y_data=1+np.array(data[:,0]),\
    np.array(data[:,1])
    Lt=np.array(x_data)*np.array(y_data)
    plt.plot(np.log(x_data),Lt,marks[j],label=labls[j])
plt.ylim(0,1.7e7)
plt.xlabel('lnt')
plt.ylabel(r't . I(t)')
plt.title('(b)')
plt.text(4,.45e7,'t.I(t) vs lnt')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)  
plt.tight_layout()
plt.show()