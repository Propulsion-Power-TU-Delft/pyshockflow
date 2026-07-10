import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
from pyshockflow import RiemannProblem
from pyshockflow import *

markerSize = 4
lw=1.5

fluid = FluidIdeal(1.4, 287.05)

# SOLUTIONS
inputFiles = [
    '../analytical/solutions/Test1.pik',
    '../roe_friction_0/Results/Test1_NX_100/Results.pik',
    '../roe_friction_0.03/Results/Test1_NX_100/Results.pik',
    '../roe_friction_0.3/Results/Test1_NX_100/Results.pik',
    '../roe_friction_1.0/Results/Test1_NX_100/Results.pik',
    '../roe_friction_3.0/Results/Test1_NX_100/Results.pik',
    ]

labels = [
    'Analytical ' + r'$(C_{\rm f} = 0)$',
    r'$C_{\rm f} = 0$',
    r'$C_{\rm f} = 0.03$',
    r'$C_{\rm f} = 0.3$',
    r'$C_{\rm f} = 1.0$',
    r'$C_{\rm f} = 3.0$',
]

# Generate colors
cmap = plt.cm.viridis
colors = cmap(np.linspace(0, 1, len(inputFiles)-1))

figSize = (10, 3.5)
fig, ax = plt.subplots(1,3, figsize=figSize)
for iInput in range(len(inputFiles)):

    with open(inputFiles[iInput], 'rb') as file:
        res = pickle.load(file)

    # Analytical solution: black reference line
    if iInput == 0:
        color = 'k'
        ax[0].plot(res.x+0.5, res.rho[:,-1], label=labels[iInput], lw=lw, color=color)
        ax[1].plot(res.x+0.5, res.u[:,-1], lw=lw, color=color)
        ax[2].plot(res.x+0.5, res.p[:,-1], lw=lw, color=color)

    else:
        color = colors[iInput-1]

        ax[0].plot(
            res['X Coords'][1:-1],
            res['Primitive']['Density'][1:-1,-1],
            label=labels[iInput],
            lw=lw,
            color=color
        )

        ax[1].plot(
            res['X Coords'][1:-1],
            res['Primitive']['Velocity'][1:-1,-1],
            lw=lw,
            color=color
        )

        ax[2].plot(
            res['X Coords'][1:-1],
            res['Primitive']['Pressure'][1:-1,-1],
            lw=lw,
            color=color
        )

ax[0].set_ylabel(r'$\rho$')
ax[1].set_ylabel(r'$u$')
ax[2].set_ylabel(r'$p$')

for axx in ax:
    axx.set_xlabel(r'$x$')
    axx.grid(alpha=.3)

fig.legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, +1.2))

plt.tight_layout()
fig.savefig('Pictures/Test1_res1.pdf', bbox_inches='tight')






figSize = (10, 3.5)
fig, ax = plt.subplots(1,3, figsize=figSize)
for iInput in range(len(inputFiles)):

    with open(inputFiles[iInput], 'rb') as file:
        res = pickle.load(file)

    # Analytical solution: black reference line
    if iInput == 0:
        color = 'k'
        
        entropy = fluid.computeEntropy_p_rho(res.p[:,-1], res.rho[:,-1])
        mach = fluid.computeMach_u_p_rho(res.u[:,-1], res.p[:,-1], res.rho[:,-1])
        total_energy = fluid.computeStaticEnergy_p_rho(res.p[:,-1], res.rho[:,-1]) + 0.5*res.u[:,-1]**2
        
        ax[0].plot(res.x+0.5, total_energy, label=labels[iInput], lw=lw, color=color)
        ax[1].plot(res.x+0.5, mach, lw=lw, color=color)
        ax[2].plot(res.x+0.5, entropy, lw=lw, color=color)

    else:
        color = colors[iInput-1]
        
        entropy = fluid.computeEntropy_p_rho(res['Primitive']['Pressure'][1:-1,-1], res['Primitive']['Density'][1:-1,-1])
        mach = fluid.computeMach_u_p_rho(res['Primitive']['Velocity'][1:-1,-1], res['Primitive']['Pressure'][1:-1,-1], res['Primitive']['Density'][1:-1,-1])
        total_energy = fluid.computeStaticEnergy_p_rho(res['Primitive']['Pressure'][1:-1,-1], res['Primitive']['Density'][1:-1,-1]) + 0.5*res['Primitive']['Velocity'][1:-1,-1]**2
        
        ax[0].plot(
            res['X Coords'][1:-1],
            total_energy,
            label=labels[iInput],
            lw=lw,
            color=color
        )

        ax[1].plot(
            res['X Coords'][1:-1],
            mach,
            lw=lw,
            color=color
        )

        ax[2].plot(
            res['X Coords'][1:-1],
            entropy,
            lw=lw,
            color=color
        )

ax[0].set_ylabel(r'$e_{\rm t}$')
ax[1].set_ylabel(r'$M$')
ax[2].set_ylabel(r'$s$')

for axx in ax:
    axx.set_xlabel(r'$x$')
    axx.grid(alpha=.3)

fig.legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, +1.2))

plt.tight_layout()
fig.savefig('Pictures/Test1_res2.pdf', bbox_inches='tight')

    
    

plt.show()
