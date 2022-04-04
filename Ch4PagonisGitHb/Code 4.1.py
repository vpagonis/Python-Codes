# The nearest neighbors distribution
# Fig 1 in Pagonis and Kulp paper
from scipy.spatial.distance import cdist 
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ne,  nh,  side,   a=\
100, 500, 100e-9, 0.11e-9
fig = plt.figure()
ax = fig.add_subplot(121, projection='3d',xlabel='x [nm]',\
ylabel='y [nm]')
pose =side* np.random.rand(ne, 3)
posh =side* np.random.rand(nh, 3)
x,y,z=1e9*pose[:,0],1e9*pose[:,1],1e9*pose[:,2]
ax.scatter(x, y, z, c='red', marker='o')
plt.title("(a)")
x,y,z=1e9*posh[:,0],1e9*posh[:,1],1e9*posh[:,2]
ax.scatter(x, y, z, c='blue', marker='^')
ax = fig.add_subplot(122)
distances=cdist(pose, posh)
rho=nh/(side**3.0)
rhoprime=(4*np.pi*rho/3)*(a**3.0)
mindistances=1e9*np.array([min(xi) for xi in distances])
plt.hist(mindistances,8)
plt.title("(b)")
plt.tight_layout() 
x=np.arange(0.0,17.0e-9,.3e-9)
rprime=(4*np.pi*rho/3)**(1/3) * x
plt.plot(x*1e9,22/0.4*(rprime**2.0)*np.exp((-rprime**3.0)),\
linewidth=2)
plt.xlabel('Nearest neighbor distance [nm]')
plt.ylabel('Distribution')
plt.tight_layout() 
plt.show()