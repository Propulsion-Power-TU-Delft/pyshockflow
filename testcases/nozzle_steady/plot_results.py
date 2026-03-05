import numpy as np
import matplotlib.pyplot as plt
import pickle
from pyshockflow.plot_styles import *
from scipy.optimize import fsolve

pressureList = [45, 75, 90, 94, 97]
pickleList = ['Results/outletPressure_%ikPa_NX_200/Results.pik' %pressure for pressure in pressureList]


fig, axes = plt.subplots(1, 2, figsize=(9, 4))

for i, pickleFile in enumerate(pickleList):
    with open(pickleFile, 'rb') as file:
        solution = pickle.load(file)
    
    if i == 0: # compute also reference
        xArea = solution['X Coords']
        areaRatio =  solution['Area'] / np.min(solution['Area'])
        gammaFluid = solution['Fluid'].gmma
        def machFunction(machLocal, areaRatioLocal, gammaFluid):
            residual = areaRatioLocal - 1/machLocal * (2/(gammaFluid+1) * (1 + (gammaFluid-1)/2 * machLocal**2))**((gammaFluid+1)/(2*(gammaFluid-1)))
            return residual

        machTheory = np.zeros(len(xArea))
        idThroat = np.argmin(areaRatio)
        for iPoint in range(len(xArea)):
            if iPoint < idThroat:
                machTheory[iPoint] = fsolve(machFunction, 0.1, args=(areaRatio[iPoint], gammaFluid))[0]
            else:
                machTheory[iPoint] = fsolve(machFunction, 1.2, args=(areaRatio[iPoint], gammaFluid))[0]

        axes[0].plot(xArea[::10], machTheory[::10], 'ko', mfc='none' ,label=r'Supersonic reference')
        
    xCoords = solution['X Coords'][1:-1]
    density = solution['Primitive']["Density"][1:-1,-1]
    pressure = solution['Primitive']["Pressure"][1:-1,-1]
    velocity = solution['Primitive']["Velocity"][1:-1,-1]
    mach = solution['Fluid'].computeMach_u_p_rho(velocity, pressure, density)
    entropy = solution['Fluid'].computeEntropy_p_rho(pressure, density)
    totalPressure = solution['Fluid'].computeTotalPressure_p_M(pressure, mach)
    temperature = solution['Fluid'].computeTemperature_p_rho(pressure, density)
    totalTemperature = solution['Fluid'].computeTotalTemperature_T_M(temperature, mach)
    
    
    axes[0].plot(xCoords, mach, label=r'$p_{\rm out}=%i$ kPa' %(pressureList[i]))
    axes[0].set_ylabel(r'$M$')
    
    axes[1].plot(xCoords, pressure/1e3)
    axes[1].set_ylabel(r'$p$ [kPa]')
    
    for ax in axes:
        ax.set_xlabel(r'$x$')
        ax.grid(alpha=.3)

fig.legend(loc='upper center', bbox_to_anchor=(0.55, 1.18), ncol=3)
plt.tight_layout()
plt.savefig('Pictures/mach_pressure_ideal_nozzle.pdf', bbox_inches='tight')
    
plt.show()
        
        
        
    