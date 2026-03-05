import numpy as np
import matplotlib.pyplot as plt
import pickle
from pyshockflow.plot_styles import *

points = [500] 
inputFiles = ['Results/rocket_NX_%i/Results.pik' %(_) for _ in points]
figsize=(4.5,3.5)

# # OPENING TRANSIENT PRESSURE DISTRIBUTION
plt.figure(figsize=figsize)
for ii, inputFile in enumerate(inputFiles):
    with open(inputFile, 'rb') as file:
        result = pickle.load(file)
    
    time = result['Time']
    nTimes = len(time)*0.00377/time[-1]
    iTimes = np.linspace(0, nTimes-1, 10, dtype=int)
    
    for it in range(len(iTimes)):
        plt.plot(result['X Coords'][1:-2], result['Primitive']['Pressure'][1:-2,iTimes[it]]/1e5, label='t/T=%.2f' %(iTimes[it]/iTimes[-1]))
plt.legend()
plt.xlabel(r'$x \ \rm{[m]}$')
plt.ylabel(r'$p \ \rm{[bar]}$')
plt.grid(alpha=.3)

# PRESSURE SENSORS IN TIME
plt.figure(figsize=figsize)
for ii, inputFile in enumerate(inputFiles):
    with open(inputFile, 'rb') as file:
        result = pickle.load(file)
    locations = [1/4, 1/2, 3/4]
    length = result['X Coords'][-1] - result['X Coords'][0]
    for location in locations:
        idx = np.argmin(np.abs(result['X Coords']-location*length))
        plt.plot(result['Time']*1e3, result['Primitive']['Pressure'][idx,:]/1e5, label='x/L=%.2f' %(location))
plt.legend()
plt.xlabel(r'$t \ \rm{[ms]}$')
plt.ylabel(r'$p \ \rm{[bar]}$')
plt.grid(alpha=.3)
plt.savefig('Pictures/sensor_pressure_rocket.pdf', bbox_inches='tight')

# MACH NUMBER AT NOZZLE THROAT AND EXIT IN TIME
plt.figure(figsize=figsize)
for ii, inputFile in enumerate(inputFiles):
    with open(inputFile, 'rb') as file:
        result = pickle.load(file)
    idxThroat = np.argmin(result['Area'])
    mach = result['Fluid'].computeMach_u_p_rho(result['Primitive']['Velocity'], result['Primitive']['Pressure'], result['Primitive']['Density'])
    plt.plot(result['Time']*1e3, mach[idxThroat,:], label='Throat')
    plt.plot(result['Time']*1e3, mach[-2,:], label='Exit')
plt.legend()
plt.xlabel(r'$t \ \rm{[ms]}$')
plt.ylabel(r'$M$')
plt.grid(alpha=.3)
plt.savefig('Pictures/nozzle_mach_rocket.pdf', bbox_inches='tight')

# MASS FLOW RATE EXIT
plt.figure(figsize=figsize)
for ii, inputFile in enumerate(inputFiles):
    with open(inputFile, 'rb') as file:
        result = pickle.load(file)
    idxThroat = np.argmin(result['Area'])
    
    # fix the initial time istant
    
    
    massflowExit = result['Primitive']['Density'][-1,:]*result['Primitive']['Velocity'][-1,:]*result['Area'][-1]
    massflowThroat = result['Primitive']['Density'][idxThroat,:]*result['Primitive']['Velocity'][idxThroat,:]*result['Area'][idxThroat]
    plt.plot(result['Time']*1e3, massflowThroat, label='Throat')
    plt.plot(result['Time']*1e3, massflowExit, label='Exit')
plt.legend()
plt.xlabel(r'$t \ \rm{[ms]}$')
plt.ylabel(r'$\dot{m} \ \rm{[kg/s]}$')
plt.grid(alpha=.3)
plt.savefig('Pictures/massflow_rocket.pdf', bbox_inches='tight')

# THRUST EXIT
plt.figure(figsize=figsize)
for ii, inputFile in enumerate(inputFiles):
    with open(inputFile, 'rb') as file:
        result = pickle.load(file)
    idxThroat = np.argmin(result['Area'])
    reactionThrust = result['Primitive']['Density'][-1,:]*(result['Primitive']['Velocity'][-1,:]**2)*result['Area'][-1]
    pressureThrust = result['Primitive']['Pressure'][-1,:]*result['Area'][-1]
    plt.plot(result['Time']*1e3, reactionThrust/1e3, label='Reaction thrust')
    plt.plot(result['Time']*1e3, pressureThrust/1e3, label='Pressure thrust')
    # plt.plot(result['Time']*1e3, (pressureThrust+reactionThrust)/1e3, label='total thrust')
plt.xlabel(r'$t \ \rm{[ms]}$')
plt.ylabel(r'$T \ \rm{[kN]}$')
plt.legend()
plt.grid(alpha=.3)
plt.savefig('Pictures/thrust_rocket.pdf', bbox_inches='tight')






plt.show()
