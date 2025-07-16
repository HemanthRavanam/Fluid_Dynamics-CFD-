import sys
import numpy as np
import scipy.linalg as sc 
import matplotlib.pyplot as plt
import math as mt
import matplotlib.pylab as pylab
import os
from chebdif import chebdif

# execute the script that in turn runs pxyst
out = os.system("bash read_stat.sh xy.stat stat_vel.in stat_vel.data")
out = os.system("bash read_stat.sh xy.stat stat_bud.in stat_bud.data")

#import subprocess
#subprocess.run(["bash", "read_stat.sh", "xy.stat", "stat_vel.in", "stat_vel.data"], check=True)


stat = np.loadtxt('stat_vel.data')
stat1 = np.loadtxt('stat_bud.data')
re=4200
y=stat[:,0]
ny=stat.shape[0]
y_,D_=chebdif(ny,1)
D=D_[0,:,:]

dudy = D@stat[:,1]
utau = np.sqrt((dudy[-1]-dudy[0])/2/re)
lstar = 1/utau/re
retau = 1/lstar
print('Re=%.5f Retau=%.5f' % (re,retau) )
yplus = y/lstar
u = stat[:,2]
labels = ['Convection', 'Production', 'Dissipation', 'Turb Diff', 'Velp', 'Visc Diff']
data = stat1[:, 1:7]  

# Non-dimensionalize each column
non_dimensional = [(col * lstar) / utau**3 for col in data.T]  

# Plotting
plt.figure(figsize=(10,6))
for i in range(6):
    plt.plot(yplus, non_dimensional[i], label=labels[i])

plt.xlabel(r'$y^+$')
plt.ylabel(r'Loss                          Gain')
plt.title('Inner scaled Budget Terms in Wall Units')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()