import sys
import numpy as np
import scipy.linalg as scla
import matplotlib.pyplot as plt
import math as mt
import matplotlib.pylab as pylab
from chebdif import chebdif

# Load data
amp = np.loadtxt('amp.stat')
re = 4200
reb = amp[0,12]
t = amp[913:,0]
n = t.size

# Compute utau and Retau
utau1 = np.sqrt(amp[913:,13] / re)
utau2 = np.sqrt(-amp[913:,14] / re)
retau1 = utau1 * re
retau2 = utau2 * re
retau = (np.mean(retau1[913:n]) + np.mean(retau2[913:n])) / 2

# Compute tau
tau1 = amp[913:,13] / re
tau2 = np.abs(amp[913:,14] / re)
tau = (np.mean(tau1[913:n]) + np.mean(tau2[913:n])) / 2

# Print 
print('Re=%.5f Reb=%.5f Retau=%.5f tau=%.5f' % (re, reb, retau, tau))

# Plot Retau
plt.figure()
plt.plot(t, retau1, label='upper wall')
plt.plot(t, retau2, label='lower wall')
plt.plot((t[0], t[-1]), (retau, retau), label='average')
plt.xlabel('time $t$')
plt.ylabel(r'$Re_{\tau}$')
plt.legend(loc='lower right')
plt.title(r'Averaged $Re_{\tau}$')

# Plot tau
plt.figure()
plt.plot(t, tau1, label='upper wall')
plt.plot(t, tau2, label='lower wall')
plt.plot((t[0], t[-1]), (tau, tau), label='average')
plt.xlabel('time $t$')
plt.ylabel(r'$\tau$')
plt.legend(loc='lower right')
plt.title(r'Averaged wall shear stress $\tau$')

# Show both figures
plt.show()