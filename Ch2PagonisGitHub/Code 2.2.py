# Solution of the OTOR using the solver solve_ivp
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
# Numerical parameters for OTOR.
N,    An,   Am,     s,   E , hr ,  kB= \
1e10, 1e-8, 1e-8,  1e12, 1 , 1,  8.617e-5
# Define derivative function
def f(t, y):
    dydt=[
		    -y[0]*s*np.exp(-E/(kB*(273+hr*t)))+y[1]*An*(N-y[0]),
		    y[0]*s*np.exp(-E/(kB*(273+\
    hr*t)))- y[1]*An*(N-y[0])-y[2]*Am*y[1],
    -y[2]*Am*y[1]]
    return dydt
tspan = np.linspace(0, 180,30)
yinit = [1e10,0,1e10]
# Solve differential equation
sol = solve_ivp(lambda t, y: f(t, y), [tspan[0], tspan[-1]],\
yinit, t_eval=tspan)
temps=hr*sol.t.reshape(-1,1)
nvals,ncvals,mvals=np.squeeze(sol.y)
# Plots
plt.plot(temps,nvals,'o-',label='n(T)')
plt.ylabel('n(T)  and   TL')
plt.ylim(0,1.2e10)
plt.title('OTOR model, using ODE solver   "solve_ivp"')
plt.xlabel(r'Temperature T [$^{o}$C]')
plt.plot(sol.t, 40*Am*ncvals*mvals,'o-',label='TLx40')
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.tight_layout()
plt.show()