# Transform CW-OSL to pseudo-LMOSL for KST4 feldspar
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
from scipy import optimize
from prettytable import PrettyTable
files=('OSL0s.txt','OSL100s.txt','OSL1000s.txt')
sym=('o-','+-','^-')
labs=('0 s of IR','100 s of IR','1000 s of IR')
plt.subplot(2,1, 1)
for i in range(0,3,1):
    data = np.loadtxt(files[i])
    x_data, y_data =0.1+np.array(data[:,0]),\
    0.1+np.array(data[:,1])
    plt.plot(np.log(x_data),y_data,sym[i],label=labs[i])
plt.ylabel('CW-OSL [cts/s]')
plt.xlabel('Time [s]')
plt.title('(a)   CW-OSL: KST4 feldspar')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.subplot(2,1,2)
for i in range(0,3,1):
    data = np.loadtxt(files[i])
    x_data, y_data=np.array(data[:,0]),\
    np.array(data[:,1])
    x_data =np.sqrt(x_data*2*len(x_data))
    y_data=x_data*y_data
    plt.plot(x_data,y_data,sym[i],label=labs[i])
plt.ylabel('I(u) [a.u.]')
plt.xlabel('u [s]')
plt.title('(b)   Transformed CW-OSL into pseudo-LM-OSL')
leg = plt.legend(loc='right')
leg.get_frame().set_linewidth(0.0)
plt.tight_layout()
plt.show()
# figure = plt.gcf() # get current figure
# # figure.set_size_inches(8, 6)
# # when saving, specify the DPI
# import os
# os.chdir('C:/Users/Bill Pagonis/Desktop')
# plt.savefig("vpplot.png", dpi = 600)
