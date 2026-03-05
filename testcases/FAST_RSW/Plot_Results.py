from pyshockflow import Driver
import numpy as np
import matplotlib.pyplot as plt
import pickle
# from PyShockflow.styles import *
from scipy.optimize import fsolve
from scipy.optimize import brentq

# INPUT
inputPkls = [
                'Results/test_coolprop_NX_250/Results.pik',
                'Results/test_coolprop_NX_500/Results.pik',
                'Results/test_coolprop_NX_1000/Results.pik',
            ]
labels = [
            '250', 
            '500',
            '1000',
            ]

linestyles = ['--', '-.', ':']
lw=2

# analytical solution
soundSpeed = 34.21  # sound speed in the left state
rswSpeed = 35       # expected speed of the rarefaction shock wave
rhoL, rhoR = 188.19, 127.27
uL, uR = 0.0, 16.78
pL, pR = 9.122E5, 8.017E5


datas = []
for inputPkl in inputPkls:    
    with open(inputPkl, 'rb') as file:
        data = pickle.load(file)
    datas.append(data)
    
# compute analytical solution at final time
time = datas[0]['Time']
lasttime = datas[0]['Time'][-1]
x_rsw = 1 - rswSpeed*lasttime
print(f"The rarefaction wave reached the exit at x = {x_rsw} [m]")

x_analytical = data['X Coords']
rho_analytical = np.zeros_like(data['Primitive']['Density'][:, -1])
u_analytical = np.zeros_like(data['Primitive']['Density'][:, -1])
p_analytical = np.zeros_like(data['Primitive']['Density'][:, -1])
for i in range(len(x_analytical)):
    x = x_analytical[i]
    if x<=x_rsw:
        rho_analytical[i] = rhoL
        u_analytical[i] = uL
        p_analytical[i] = pL
    else:
        rho_analytical[i] = rhoR
        u_analytical[i] = uR
        p_analytical[i] = pR


fig, axes = plt.subplots(1, 3, figsize=(10, 3.5))
(ax1, ax2, ax3) = axes

# --- Density ---
ax1.plot(x_analytical, rho_analytical, linestyle='-', lw=lw, label='Reference')
for data, label, ls in zip(datas, labels, linestyles):
    ax1.plot(data['X Coords'], data['Primitive']['Density'][:, -1], ls, label=label)
ax1.set_xlabel(r'$x$ [m]')
ax1.set_ylabel(r'$\rho$ [kg/m$^3$]')
fig.legend(loc='upper center', ncol=len(labels)+2, bbox_to_anchor=(0.53, 1.1))

# --- Velocity (middle) ---
ax2.plot(x_analytical, u_analytical, linestyle='-', lw=lw)
for data, label, ls in zip(datas, labels, linestyles):
    ax2.plot(data['X Coords'], data['Primitive']['Velocity'][:, -1], ls, label=label)
ax2.set_xlabel(r'$x$ [m]')
ax2.set_ylabel(r'$u$ [m/s]')

# --- Pressure (last) ---
ax3.plot(x_analytical, p_analytical/1E5, linestyle='-', lw=lw)
for data, label, ls in zip(datas, labels, linestyles):
    ax3.plot(data['X Coords'], data['Primitive']['Pressure'][:, -1] / 1e5, ls, label=label)
ax3.set_xlabel(r'$x$ [m]')
ax3.set_ylabel(r'$p$ [bar]')

for ax in axes:
    ax.grid(alpha=0.2)
    ax.set_xlim(0.25, 0.7)

fig.tight_layout()
plt.savefig('FAST_RSW_Allresults.pdf', bbox_inches='tight')



plt.show()