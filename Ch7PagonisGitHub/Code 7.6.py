# Peak shape transformation of ITL signals
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from prettytable import PrettyTable
fils=['LBOITL1.txt','LBOITL3.txt','LBOITL6.txt',
'LBOITL7.txt','LBOITL8.txt'] 
plt.subplot(1,2,1)
marks=['+','o','^','.','x']
labls=[r'95$^o$C',r'120$^o$C',r'145$^o$C',r'170$^o$C',\
r'195$^o$C']
for j in range(2,5):
    data = np.loadtxt(fils[j],delimiter=',')
    x_data,y_data =1+np.array(data[:,0]), np.array(data[:,1])
    plt.plot(x_data,y_data,marks[j]+'-',label=labls[j])
    plt.xlabel('Stimulation time t [s]')
    plt.ylabel(r'ITL signal  [a.u.]')
plt.title('(a)')
plt.text(300,.7,'ITL signal')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0) 
plt.subplot(1,2,2)
for j in range(2,5):
    data = np.loadtxt(fils[j],delimiter=',')
    x_data, y_data = 1+np.array(data[:, 0]), np.array(data[:, 1])
    Lt=np.array(x_data)*np.array(y_data)
    plt.plot(np.log(x_data),Lt,marks[j],label=labls[j])
    plt.xlabel('lnt')
    plt.ylabel(r't . I(t)')
plt.ylim(0,60)
plt.title('(b)')
plt.text(2.2,42,'t.I(t) vs lnt')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)  
plt.tight_layout()
plt.show()