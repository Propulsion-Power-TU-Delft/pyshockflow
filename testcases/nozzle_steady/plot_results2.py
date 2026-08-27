import numpy as np
import matplotlib.pyplot as plt
import pickle
from pyshockflow.fluid import FluidIdeal
from pyshockflow.plot_styles import *
from scipy.optimize import fsolve

pressureList = [45, 75, 90, 95, 97]
colors = plt.cm.viridis(np.linspace(1, 0, len(pressureList)))
pickleList = ['Results/outletPressure_%ikPa_NX_250/Results.pik' %pressure for pressure in pressureList]
fluid = FluidIdeal(1.4, 287.05)

fig, axes = plt.subplots(1, 2, figsize=(9, 4))

for i, pickleFile in enumerate(pickleList):
    with open(pickleFile, 'rb') as file:
        solution = pickle.load(file)
        
    xCoords = solution['X Coords'][1:-1]
    density = solution['Primitive']["Density"][1:-1,-1]
    pressure = solution['Primitive']["Pressure"][1:-1,-1]
    velocity = solution['Primitive']["Velocity"][1:-1,-1]
    mach = fluid.computeMach_u_p_rho(velocity, pressure, density)
    entropy = fluid.computeEntropy_p_rho(pressure, density)
    totalPressure = fluid.computeTotalPressure_p_M(pressure, mach)
    temperature = fluid.computeTemperature_p_rho(pressure, density)
    totalTemperature = fluid.computeTotalTemperature_T_M(temperature, mach)
    
    skip = 1
    axes[0].plot(xCoords[::skip], mach[::skip], '--', color=colors[i], mfc = 'none', ms=5)
    axes[0].set_ylabel(r'$M$')
    
    axes[1].plot(xCoords[::skip], pressure[::skip]/1e3, '--', color=colors[i], mfc = 'none', ms=5)
    axes[1].set_ylabel(r'$p$ [kPa]')
    
    for ax in axes:
        ax.set_xlabel(r'$x$')
        ax.grid(alpha=.3)
    
analyticalFile = 'nozzle_solutions_analytical.pkl'
with open(analyticalFile, 'rb') as file:
    analyticalSolutions = pickle.load(file)

for i, p_out in enumerate(pressureList):
    sol = analyticalSolutions[p_out*1e3]
    xCoords = sol['x']
    mach = sol['M']
    pressure = sol['p']
    
    skip = 1
    axes[0].plot(xCoords[::skip], mach[::skip], '-', color=colors[i], mfc = 'none', label=r'$p_{\rm out}=%i$ kPa' %(pressureList[i]))
    axes[1].plot(xCoords[::skip], pressure[::skip]/1e3, '-', color=colors[i], mfc = 'none')

fig.legend(loc='upper center', bbox_to_anchor=(0.55, 1.18), ncol=3)
plt.tight_layout()
plt.savefig('Pictures/mach_pressure_ideal_nozzle.pdf', bbox_inches='tight')
    
plt.show()
        
        
        
    