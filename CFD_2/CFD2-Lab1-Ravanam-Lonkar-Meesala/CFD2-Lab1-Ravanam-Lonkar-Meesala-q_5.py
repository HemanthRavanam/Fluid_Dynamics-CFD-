import sys
import numpy as np
import scipy.linalg as scla
import matplotlib.pyplot as plt
import math as mt
import matplotlib.pylab as pylab
from chebdif import chebdif

amp = np.loadtxt('amp2.stat')
re=4200
reb=amp[0,12]

t=amp[:,0]
n=t.size
utau1 = np.sqrt( amp[:,13]/re)
utau2 = np.sqrt(-amp[:,14]/re)
retau1 = utau1*re
retau2 = utau2*re

retau = (np.mean(retau1[int(n/2):n])+np.mean(retau2[int(n/2):n]))/2
print('Re=%.5f Reb=%.5f Retau=%.5f' % (re,reb,retau) )

mask = t > 200
t_masked = t[mask]

utau1_masked = utau1[mask]
utau2_masked = utau2[mask]

tauwall1 = np.square(utau1_masked)
tauwall2 = np.square(utau2_masked)
tauwall_averaged = (tauwall1 + tauwall2) / 2
 
q_dash = np.gradient(tauwall_averaged, t_masked)      
denominator = np.mean(q_dash**2)
    
N = len(q_dash)                           
s = 200
rho = np.zeros(s)                         

for i in range(s):                      
    n_overlap = N - i                      
    numerator = np.sum(q_dash[:n_overlap] * q_dash[i:])  
    rho[i] = numerator / (denominator * N)  

lags = np.arange(s) 

plt.figure(figsize=(8, 4))
plt.plot(lags, rho, linestyle='-', color='blue', label=r'$\rho(s)$')
plt.xlabel('s')
plt.ylabel(r'$\rho(s)$')
plt.title(r'$\rho(s)$ vs. Lag $s$')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

