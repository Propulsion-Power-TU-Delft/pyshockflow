import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
from pyshockflow import RiemannProblem
from pyshockflow import *
from pyshockflow.fluid import FluidIdeal

markerSize = 4
lw=1.5

fluid = FluidIdeal(1.4, 287.05)
R = fluid.Rgas
gmma = fluid.gmma
cp = (gmma*R)/(gmma-1)

# SOLUTIONS
inputFiles = [
    '../heatflux_00000/Results/heatflux_NX_300/Results.pik',
    '../heatflux_01000/Results/heatflux_NX_300/Results.pik',
    '../heatflux_05000/Results/heatflux_NX_300/Results.pik',
    '../heatflux_10000/Results/heatflux_NX_300/Results.pik',
    ]

labels = [
    r'$\dot{q}_{\rm w} = 0~\rm{kW/m^2}$',
    r'$\dot{q}_{\rm w} = 1~\rm{kW/m^2}$',
    r'$\dot{q}_{\rm w} = 5~\rm{kW/m^2}$',
    r'$\dot{q}_{\rm w} = 10~\rm{kW/m^2}$',
    ]

heatFlux = [
    0.000,
    1000.000,
    5000.000,
    10000.000
]

# Generate colors
cmap = plt.cm.viridis
colors = cmap(np.linspace(0, 1, len(inputFiles)))

figSize = (8, 3.5)
fig, ax = plt.subplots(1,2, figsize=figSize)
for iInput in range(len(inputFiles)):

    with open(inputFiles[iInput], 'rb') as file:
        res = pickle.load(file)
    
    mach = fluid.computeMach_u_p_rho(res['Primitive']['Velocity'][1:-1,-1], res['Primitive']['Pressure'][1:-1,-1], res['Primitive']['Density'][1:-1,-1])
    rho = res['Primitive']['Density'][1:-1,-1]
    u = res['Primitive']['Velocity'][1:-1,-1]
    T = fluid.computeTemperature_p_rho(res['Primitive']['Pressure'][1:-1,-1], res['Primitive']['Density'][1:-1,-1])
    totalT = T * (1 + (gmma-1)/2 * mach**2)
    area = 0.01
    length = 10.0
    mach_exit = mach[-1]
    R = np.sqrt(area/np.pi)
    p_outlet = 90000
    p = res['Primitive']['Pressure'][1:-1,-1]
    totalP = p * (1 + (fluid.gmma-1)/2 * mach**2)**(fluid.gmma/(fluid.gmma-1))

    color = colors[iInput]

    ax[0].plot(
        res['X Coords'][1:-1]/length,
        totalT/totalT[0],
        '--',
        lw=lw,
        color=color
    )

    ax[1].plot(
        res['X Coords'][1:-1]/length,
        mach,
        '--',
        lw=lw,
        color=color
    )
    
    # global check of integral quantites
    deltaTt = totalT[-1] - totalT[0]
    totalQ = heatFlux[iInput] * 2*np.pi*R * length
    mdot = rho[5]*u[5]*area
    expectedDeltaT = totalQ / (mdot * cp)
    print(f"Delta Tt for qw={heatFlux[iInput]}: {deltaTt:.4f} K, expected value: {expectedDeltaT:.4f} K")
    
    deltaPt = totalP[-1] - totalP[0]
    print(f"Delta Pt for qw={heatFlux[iInput]}: {deltaPt/totalP[0]:.4f} %")


# analytical solutions
analyticalFile = '../analytical/analytical_results.pkl'
with open(analyticalFile, 'rb') as file:
    analyticalSolutions = pickle.load(file)

for iInput in range(len(inputFiles)):
    sol = analyticalSolutions[heatFlux[iInput]]
    qw = sol['qw']
    xCoords = sol['x']
    mach = sol['M']
    
    color = colors[iInput]
    
    ax[0].plot(
        sol['x']/length,
        sol['Tt']/sol['Tt'][0],
        '-',
        color=color,
        lw=lw,
        label=labels[iInput],
    )
    
    ax[1].plot(
        sol['x']/length,
        sol['M'],
        '-',
        color=color,
        lw=lw,
    )


ax[0].set_ylabel(r'$T_{\rm t}/T_{\rm t, in}$')
ax[1].set_ylabel(r'$M$')

for axx in ax:
    axx.set_xlabel(r'$x/L$')
    axx.grid(alpha=.3)

fig.legend(loc='upper center', ncol=2, bbox_to_anchor=(0.5, +1.23))
fig.tight_layout()
fig.savefig('Pictures/validation_heatflux_duct.pdf', bbox_inches='tight')


plt.show()
