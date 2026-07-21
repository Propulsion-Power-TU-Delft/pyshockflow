import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
from pyshockflow import RiemannProblem
from pyshockflow import *

markerSize = 4
lw=1.5

fluid = FluidIdeal(1.4, 287.05)
R = fluid.Rgas
gmma = fluid.gmma
cp = (gmma*R)/(gmma-1)

# SOLUTIONS
inputFiles = [
    '../heatflux_00000/Results/heatflux_NX_100/Results.pik',
    '../heatflux_01000/Results/heatflux_NX_100/Results.pik',
    '../heatflux_05000/Results/heatflux_NX_100/Results.pik',
    '../heatflux_10000/Results/heatflux_NX_100/Results.pik',
    ]

labels = [
    r'$\dot{q}_{\rm w} = 0 \ \rm{[W/m^2]}$',
    r'$\dot{q}_{\rm w} = 1000 \ \rm{[W/m^2]}$',
    r'$\dot{q}_{\rm w} = 5000 \ \rm{[W/m^2]}$',
    r'$\dot{q}_{\rm w} = 10000 \ \rm{[W/m^2]}$',
    ]

qw = [
    0.000,
    1000.000,
    5000.000,
    10000.000
]

# Generate colors
cmap = plt.cm.viridis
colors = cmap(np.linspace(0, 1, len(inputFiles)))

figSize = (12, 3.5)
fig, ax = plt.subplots(1,3, figsize=figSize)
for iInput in range(len(inputFiles)):

    with open(inputFiles[iInput], 'rb') as file:
        res = pickle.load(file)
    
    mach = fluid.computeMach_u_p_rho(res['Primitive']['Velocity'][1:-1,-1], res['Primitive']['Pressure'][1:-1,-1], res['Primitive']['Density'][1:-1,-1])
    rho = res['Primitive']['Density'][1:-1,-1]
    u = res['Primitive']['Velocity'][1:-1,-1]
    dmdx = np.gradient(mach, res['X Coords'][1:-1])
    T = fluid.computeTemperature_p_rho(res['Primitive']['Pressure'][1:-1,-1], res['Primitive']['Density'][1:-1,-1])
    totalT = T * (1 + (gmma-1)/2 * mach**2)
    beta = mach/2 * (1.4*mach**2 + 1)/(1-mach**2)
    area = 0.01
    length = 10.0
    mach_exit = mach[-1]
    R = np.sqrt(area/np.pi)
    dmdx_analytical = beta / (rho*u*area*cp*T) * qw[iInput] * 2*np.pi*R
    skip = 8
    p_outlet = 90000
    p = res['Primitive']['Pressure'][1:-1,-1]
    totalP = p * (1 + (fluid.gmma-1)/2 * mach**2)**(fluid.gmma/(fluid.gmma-1))

    color = colors[iInput]

    ax[0].plot(
        res['X Coords'][1:-1]/length,
        totalT/totalT[0],
        label=labels[iInput],
        lw=lw,
        color=color
    )

    ax[1].plot(
        res['X Coords'][1:-1]/length,
        mach,
        lw=lw,
        color=color
    )

    ax[2].plot(
        res['X Coords'][1:-1][::skip]/length,
        dmdx_analytical[::skip]*length/mach_exit,
        'o',
        color=color,
        mfc='none'
    )
    
    ax[2].plot(
        res['X Coords'][1:-1]/length,
        dmdx*length/mach_exit,
        lw=lw,
        color=color
    )
    
    # global check of integral quantites
    deltaTt = totalT[-1] - totalT[0]
    totalQ = qw[iInput] * 2*np.pi*R * length
    mdot = rho[5]*u[5]*area
    expectedDeltaT = totalQ / (mdot * cp)
    print(f"Delta Tt for qw={qw[iInput]}: {deltaTt:.4f} K, expected value: {expectedDeltaT:.4f} K")
    
    deltaPt = totalP[-1] - totalP[0]
    print(f"Delta Pt for qw={qw[iInput]}: {deltaPt/totalP[0]:.4f} %")

ax[0].set_ylabel(r'$T_{\rm t}/T_{\rm t, in}$')
ax[1].set_ylabel(r'$M$')
ax[2].set_ylabel(r'$\frac{dM}{dx} \cdot  \frac{L}{M_{\rm out}}$')

for axx in ax:
    axx.set_xlabel(r'$x/L$')
    axx.grid(alpha=.3)

fig.legend(loc='upper center', ncol=2, bbox_to_anchor=(0.5, +1.2))
fig.tight_layout()
fig.savefig('Pictures/validation_heatflux_duct.pdf', bbox_inches='tight')


plt.show()
