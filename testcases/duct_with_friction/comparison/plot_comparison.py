import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
from pyshockflow import RiemannProblem
from pyshockflow import *

markerSize = 4
lw=1.5

fluid = FluidIdeal(1.4, 287.05)
cp = fluid.Rgas*fluid.gmma/(fluid.gmma-1)


# SOLUTIONS
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

figSize = (12, 3.5)
fig, ax = plt.subplots(1,3, figsize=figSize)
for iInput in range(len(inputFiles)):

    with open(inputFiles[iInput], 'rb') as file:
        res = pickle.load(file)
    
    mach = fluid.computeMach_u_p_rho(res['Primitive']['Velocity'][:,-1], res['Primitive']['Pressure'][:,-1], res['Primitive']['Density'][:,-1])
    dmdx = np.gradient(mach, res['X Coords'][:])
    alpha = 1.4 * mach**3 * (1 + 0.2*mach**2) / (1 - mach**2) # analytical from D'Agostino slides for Fanno flows
    area = 0.01
    length = 10.0
    mach_exit = mach[-1]
    dmdx_analytical = alpha * cf[iInput] / np.sqrt(area/np.pi)
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
        label=labels[iInput],
        lw=lw,
        color=color
    )

    ax[1].plot(
        res['X Coords'][:]/length,
        mach,
        lw=lw,
        color=color
    )

    ax[2].plot(
        res['X Coords'][:][::skip]/length,
        dmdx_analytical[::skip]*length/mach_exit,
        'o',
        color=color,
        mfc='none'
    )
    
    ax[2].plot(
        res['X Coords'][:]/length,
        dmdx*length/mach_exit,
        lw=lw,
        color=color
    )
    
    # global check of integral quantites
    deltaTt = totalT[-1] - totalT[0]
    print(f"Delta Tt for Cf={cf[iInput]}: {deltaTt/totalT.mean()*100:.4g} [%], expected value: {0:.4g} [%]")
    
    deltaPt = totalP[-1] - totalP[0]
    print(f"Delta Pt for Cf={cf[iInput]}: {deltaPt/totalP.mean()*100:.4g} [%]")
print()

ax[0].set_ylabel(r'$p_{\rm t}/p_{\rm t, in}$')
ax[1].set_ylabel(r'$M$')
ax[2].set_ylabel(r'$\frac{dM}{dx} \cdot  \frac{L}{M_{\rm out}}$')

for axx in ax:
    axx.set_xlabel(r'$x/L$')
    axx.grid(alpha=.3)

fig.legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, +1.2))
fig.tight_layout()
fig.savefig('Pictures/validation_friction_duct.pdf', bbox_inches='tight')













figSize = (12, 3.5)
fig, ax = plt.subplots(1,3, figsize=figSize)
for iInput in range(len(inputFiles)):

    with open(inputFiles[iInput], 'rb') as file:
        res = pickle.load(file)
    
    mach = fluid.computeMach_u_p_rho(res['Primitive']['Velocity'][:,-1], res['Primitive']['Pressure'][:,-1], res['Primitive']['Density'][:,-1])
    dmdx = np.gradient(mach, res['X Coords'][:])
    alpha = 1.4 * mach**3 * (1 + 0.2*mach**2) / (1 - mach**2) # analytical from D'Agostino slides for Fanno flows
    area = 0.01
    length = 10.0
    mach_exit = mach[-1]
    dmdx_analytical = alpha * cf[iInput] / np.sqrt(area/np.pi)
    skip = 8
    p = res['Primitive']['Pressure'][:,-1]
    totalP = p * (1 + (fluid.gmma-1)/2 * mach**2)**(fluid.gmma/(fluid.gmma-1))
    p_outlet = 90000
    T = fluid.computeTemperature_p_rho(res['Primitive']['Pressure'][:,-1], res['Primitive']['Density'][:,-1])
    totalT = T * (1 + (fluid.gmma-1)/2 * mach**2)
    enthalpy = T*cp
    kin_energy = 0.5*res['Primitive']['Velocity'][:,-1]**2
    enthalpy_tot = enthalpy + kin_energy

    color = colors[iInput]

    ax[0].plot(
        res['X Coords'][:]/length,
        enthalpy,
        label=labels[iInput],
        lw=lw,
        color=color
    )

    ax[1].plot(
        res['X Coords'][:]/length,
        kin_energy,
        lw=lw,
        color=color
    )
    
    ax[2].plot(
        res['X Coords'][:]/length,
        enthalpy + kin_energy,
        lw=lw,
        color=color
    )
    
    # global check of integral quantites
    deltaHt = enthalpy_tot[-1] - enthalpy_tot[0]
    print(f"Delta Ht for Cf={cf[iInput]}: {deltaHt/enthalpy_tot.mean()*100:.4g} [%], expected value: {0:.4g} [%]")

ax[0].set_ylabel(r'$h \ \rm [J/kg]$')
ax[1].set_ylabel(r'$E_{\rm kin} \ \rm [J/kg]$')
ax[2].set_ylabel(r'$h_{\rm tot} \ \rm [J/kg]$')

for axx in ax:
    axx.set_xlabel(r'$x/L$')
    axx.grid(alpha=.3)

fig.legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, +1.2))
fig.tight_layout()
fig.savefig('Pictures/validation_friction_duct_enthalpy.pdf', bbox_inches='tight')


















plt.show()
