
import sys
import numpy as np
import scipy.linalg as scla
import matplotlib.pyplot as plt
import math as mt
import matplotlib.pylab as pylab
import os
from chebdif import chebdif

# execute the script that in turn runs pxyst
out = os.system("bash read_stat.sh xy.stat stat_vel.in stat_vel.data")

stat = np.loadtxt('stat_vel.data')
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

# Load the data
stat = np.loadtxt('stat_vel.data')

# Assign each column to a named variable
y        = stat[:, 0]
u        = stat[:, 1]
v        = stat[:, 2]
w        = stat[:, 3]
u_rms    = stat[:, 4]
v_rms    = stat[:, 5]
w_rms    = stat[:, 6]
uv_dash  = stat[:, 7]
uw_dash  = stat[:, 8]
vw_dash  = stat[:, 9]

y_plus = y / lstar
u_dash = np.sqrt(u_rms)
v_dash = np.sqrt(v_rms)
w_dash = np.sqrt(w_rms)

uu_plus = u_dash * u_dash / utau**2
vv_plus = v_dash * v_dash / utau**2
ww_plus = w_dash * w_dash / utau**2 
uv_plus = uv_dash / utau**2
uw_plus = uw_dash / utau**2
vw_plus = vw_dash / utau**2
u_rms_plus = u_rms / utau
v_rms_plus = v_rms / utau    
w_rms_plus = w_rms / utau

fig, axs = plt.subplots(3, 3, figsize=(14, 10), sharex=True)
axs = axs.flatten()

# Plotting
plot_data = [
    (uu_plus, r"$\langle u'^2 \rangle^+$"),
    (vv_plus, r"$\langle v'^2 \rangle^+$"),
    (ww_plus, r"$\langle w'^2 \rangle^+$"),
    (uv_plus, r"$\langle u'v' \rangle^+$"),
    (uw_plus, r"$\langle u'w' \rangle^+$"),
    (vw_plus, r"$\langle v'w' \rangle^+$"),
    (u_rms_plus, r"$u'_{\mathrm{rms}}^+$"),
    (v_rms_plus, r"$v'_{\mathrm{rms}}^+$"),
    (w_rms_plus, r"$w'_{\mathrm{rms}}^+$"),
]

for i, (data, label) in enumerate(plot_data):
    axs[i].plot(y_plus, data, label=label)
    axs[i].set_title(label, fontsize=11)
    axs[i].grid(True, linestyle='--', linewidth=0.5)
    axs[i].legend(fontsize=9)
    axs[i].set_ylabel("Value")
    if i >= 6:
        axs[i].set_xlabel(r"$y^+$")

plt.suptitle("Turbulence Statistics in Wall Units", fontsize=14)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()

u_centreline = u[32]
u_bulk = (2 / 3) * u_centreline
uu_outer = u_dash * u_dash / u_bulk**2
vv_outer = v_dash * v_dash / u_bulk**2      
ww_outer = w_dash * w_dash / u_bulk**2
uv_outer = uv_dash / u_bulk**2
uw_outer = uw_dash / u_bulk**2
vw_outer = vw_dash / u_bulk**2
u_rms_outer = u_rms / u_bulk
v_rms_outer = v_rms / u_bulk
w_rms_outer = w_rms / u_bulk


fig, axs = plt.subplots(3, 3, figsize=(14, 10), sharex=True)
axs = axs.flatten()

# Plotting
plot_data = [
    (uu_outer, r"$\langle u'^2 \rangle$"),
    (vv_outer, r"$\langle v'^2 \rangle$"),
    (ww_outer, r"$\langle w'^2 \rangle$"),
    (uv_outer, r"$\langle u'v' \rangle$"),
    (uw_outer, r"$\langle u'w' \rangle$"),
    (vw_outer, r"$\langle v'w' \rangle$"),
    (u_rms_outer, r"$u'_{\mathrm{rms}}$"),
    (v_rms_outer, r"$v'_{\mathrm{rms}}$"),
    (w_rms_outer, r"$w'_{\mathrm{rms}}$"),
]

for i, (data, label) in enumerate(plot_data):
    axs[i].plot(y, data, label=label)
    axs[i].set_title(label, fontsize=11)
    axs[i].grid(True, linestyle='--', linewidth=0.5)
    axs[i].legend(fontsize=9)
    axs[i].set_ylabel("Value")
    if i >= 6:
        axs[i].set_xlabel(r"$y$")

plt.suptitle("Turbulence Statistics in Outer Units", fontsize=14)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()





