import matplotlib.pyplot as plt
import numpy as np
import pickle
from pyshockflow import RiemannProblem
from pyshockflow import Driver
import os
from pyshockflow.fluid import FluidIdeal

fluid = FluidIdeal(1.4, 287.05)

refinement_locations = [0.4, 0.6]

resultsFile = [
    "Results/normal_NX_500/Results.pik",
    "Results/refined_NX_500/Results.pik",
    "Results/refined_adaptive_NX_504/Results.pik"
    ]

labels = [
          'Normal', 
          'Refined', 
          'Refined Adaptive',
          ]

lss = [
    '-',
    '--',
    '-.',
]

lw = 2.5

outFolder = 'Pictures'
os.makedirs(outFolder, exist_ok=True)

resultsPickle = []
for i in range(len(resultsFile)):
    with open(resultsFile[i], 'rb') as file:
        resultsPickle.append(pickle.load(file))
        
time = resultsPickle[0]['Time']
ntimes = len(time)    
        
fig, ax = plt.subplots(2, 2, figsize=(16,10))
    
for i,results in enumerate(resultsPickle):
    results['Primitive']["Energy"] = fluid.computeStaticEnergy_p_rho(results['Primitive']["Pressure"], results['Primitive']["Density"])
    
    ax[0,0].plot(results['X Coords'], results['Primitive']["Density"][:, -1], label=labels[i], ls=lss[i], linewidth=lw)
    ax[0,1].plot(results['X Coords'], results['Primitive']["Velocity"][:, -1], label=labels[i], ls=lss[i], linewidth=lw)
    ax[1,0].plot(results['X Coords'], results['Primitive']["Pressure"][:, -1], label=labels[i], ls=lss[i], linewidth=lw)
    ax[1,1].plot(results['X Coords'], results['Primitive']["Energy"][:, -1], label=labels[i], ls=lss[i], linewidth=lw)

ax[0, 0].set_ylabel(r'Density')
ax[0, 1].set_ylabel(r'Velocity')
ax[1,  0].set_ylabel(r'Pressure')
ax[1,  1].set_ylabel(r'Energy')

# Add legend only once for the figure
fig.legend(labels, loc='upper center', ncol=3)

for row in ax:
    for col in row:
        for location in refinement_locations:
            col.axvline(location, color='r', ls=':', lw=2)
        col.set_xlabel('x')
        col.grid(alpha=.3)

fig.tight_layout()
plt.savefig(outFolder + '/Comparison.pdf', bbox_inches='tight')



plt.show()