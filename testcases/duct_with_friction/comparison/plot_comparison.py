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
cp = fluid.Rgas*fluid.gmma/(fluid.gmma-1)


# SOLUTIONS from pyshock flow
inputFiles = [
    '../roe_friction_0.000/Results/frictional_duct_NX_100/Results.pik',
    '../roe_friction_0.001/Results/frictional_duct_NX_100/Results.pik',
    '../roe_friction_0.003/Results/frictional_duct_NX_100/Results.pik',
    '../roe_friction_0.010/Results/frictional_duct_NX_100/Results.pik',
    '../roe_friction_0.030/Results/frictional_duct_NX_100/Results.pik',
    ]

labels = [
    r'$C_{\rm f} = 0.000$',
    r'$C_{\rm f} = 0.001$',
    r'$C_{\rm f} = 0.003$',
    r'$C_{\rm f} = 0.010$',
    r'$C_{\rm f} = 0.030$',
]

cf = [
    0.000,
    0.001,
    0.003,
    0.010,
    0.030,
]

# Generate colors
cmap = plt.cm.viridis
colors = cmap(np.linspace(0, 1, len(inputFiles)))

figSize = (8, 3.5)
fig, ax = plt.subplots(1,2, figsize=figSize)
for iInput in range(len(inputFiles)):

    with open(inputFiles[iInput], 'rb') as file:
        res = pickle.load(file)
    
    mach = fluid.computeMach_u_p_rho(res['Primitive']['Velocity'][:,-1], res['Primitive']['Pressure'][:,-1], res['Primitive']['Density'][:,-1])
    dmdx = np.gradient(mach, res['X Coords'][:])
    area = 0.01
    length = 10.0
    mach_exit = mach[-1]
    skip = 8
    p = res['Primitive']['Pressure'][:,-1]
    totalP = p * (1 + (fluid.gmma-1)/2 * mach**2)**(fluid.gmma/(fluid.gmma-1))
    p_outlet = 90000
    T = fluid.computeTemperature_p_rho(res['Primitive']['Pressure'][:,-1], res['Primitive']['Density'][:,-1])
    totalT = T * (1 + (fluid.gmma-1)/2 * mach**2)

    color = colors[iInput]

    ax[0].plot(
        res['X Coords'][:]/length,
        totalP/totalP[0],
        '--',
        lw=lw,
        color=color
    )

    ax[1].plot(
        res['X Coords'][:]/length,
        mach,
        '--',
        lw=lw,
        color=color
    )
    
    # global check of integral quantites
    deltaTt = totalT[-1] - totalT[0]
    print(f"Delta Tt for Cf={cf[iInput]}: {deltaTt/totalT.mean()*100:.4g} [%], expected value: {0:.4g} [%]")
    
    deltaPt = totalP[-1] - totalP[0]
    print(f"Delta Pt for Cf={cf[iInput]}: {deltaPt/totalP.mean()*100:.4g} [%]")



# analytical solutions
analyticalFile = '../analytical/results.pkl'
with open(analyticalFile, 'rb') as file:
    analyticalSolutions = pickle.load(file)

for iInput in range(len(inputFiles)):
    sol = analyticalSolutions[iInput]
    cf = sol['Cf']
    xCoords = sol['x']
    mach = sol['M']
    
    color = colors[iInput]
    
    ax[0].plot(
        sol['x']/length,
        sol['pt']/sol['pt'][0],
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
        lw=lw
    )








ax[0].set_ylabel(r'$p_{\rm t}/p_{\rm t, in}$')
ax[1].set_ylabel(r'$M$')

for axx in ax:
    axx.set_xlabel(r'$x/L$')
    axx.grid(alpha=.3)

fig.legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, +1.2))
fig.tight_layout()
fig.savefig('Pictures/validation_friction_duct.pdf', bbox_inches='tight')




















plt.show()
