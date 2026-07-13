import matplotlib.pyplot as plt
import numpy as np
import pickle
from pyshockflow import RiemannProblem
from pyshockflow import Driver
import os
from pyshockflow.fluid import FluidIdeal

fluid = FluidIdeal(1.4, 287.05)


resultsFile = ["Results/normal_NX_500/Results.pik",
               "Results/refined_NX_498/Results.pik"]

labels = ['Normal', 'Local Refinement']

outFolder = 'Pictures'
os.makedirs(outFolder, exist_ok=True)

resultsPickle = []
for i in range(len(resultsFile)):
    with open(resultsFile[i], 'rb') as file:
        resultsPickle.append(pickle.load(file))
        
fig, ax = plt.subplots(2, 2, figsize=(16,10))

for i,results in enumerate(resultsPickle):
    ax[0,0].plot(results['X Coords'][1:-1], results['Primitive']["Density"][1:-1, -1], label=labels[i])
ax[0, 0].set_ylabel(r'Density')

for i,results in enumerate(resultsPickle):
    ax[0,1].plot(results['X Coords'][1:-1], results['Primitive']["Velocity"][1:-1, -1], label=labels[i])
ax[0, 1].set_ylabel(r'Velocity')

for i,results in enumerate(resultsPickle):
    ax[1,0].plot(results['X Coords'][1:-1], results['Primitive']["Pressure"][1:-1, -1], label=labels[i])
ax[1, 0].set_ylabel(r'Pressure')



for i,results in enumerate(resultsPickle):
    results['Primitive']["Mach"] = fluid.computeMach_u_p_rho(results['Primitive']["Velocity"], results['Primitive']["Pressure"], results['Primitive']["Density"])
    ax[1,1].plot(results['X Coords'][1:-1], results['Primitive']["Mach"][1:-1, -1], label=labels[i])
ax[1, 1].set_ylabel(r'Mach')

# Add legend only once for the figure
fig.legend(labels, loc='upper center', ncol=2)

for row in ax:
        for col in row:
            col.set_xlabel('x')
            col.grid(alpha=.3)

plt.savefig(outFolder + '/Comparison.pdf', bbox_inches='tight')



plt.show()