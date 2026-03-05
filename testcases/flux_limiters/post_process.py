import matplotlib.pyplot as plt
import numpy as np
import pickle
from pyshockflow import RiemannProblem
from pyshockflow import Driver
import os


resultsFile = [
                "Results/1stOrder_NX_100/Results.pik",
                "Results/vanalbada_NX_100/Results.pik",
                "Results/minmod_NX_100/Results.pik",
                "Results/vanleer_NX_100/Results.pik",
                "Results/superbee_NX_100/Results.pik",
               ]

labels = ['Reference', # this is needed for the ref results, dont delete
          '1st order',
          'Van Albada',
          'MinMod',
          'Van Leer',
          'Superbee',
         ]

lw = 2.0
outFolder = 'Pictures'
os.makedirs(outFolder, exist_ok=True)

resultsPickle = []
for i in range(len(resultsFile)):
    with open(resultsFile[i], 'rb') as file:
        resultsPickle.append(pickle.load(file))
        
def normalize(f):
    return f
    


fig, ax = plt.subplots(1, 3, figsize=(10, 3.5))

# REFERENCE RESULTS
with open('reference.pik', 'rb') as file:
    data = pickle.load(file)
mach = data.u[:,-1]/(data.p[:,-1]*1.4/data.rho[:,-1])**0.5
ax[0].plot(data.x+0.5, data.rho[:,-1], '-k', label='Reference', lw=1.7)
# ax[0].plot(data.x+0.5, data.rho[:,-1], 'ko', mfc='none', label='Reference', ms=3.0, markevery=2)
ax[1].plot(data.x+0.5, data.u[:,-1], '-k', label='Reference', lw=1.7)
ax[2].plot(data.x+0.5, data.p[:,-1], '-k', label='Reference', lw=1.7)



colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5']
ls = ['-', '--', ':', '-.', '--', '-.']

for i,results in enumerate(resultsPickle):
    ax[0].plot(results['X Coords'][1:-1], normalize(results['Primitive']["Density"][1:-1, -1]), label=labels[i+1], color=colors[i], ls=ls[i], lw=lw)
ax[0].set_ylabel(r'$\rho$')

for i,results in enumerate(resultsPickle):
    ax[1].plot(results['X Coords'][1:-1], normalize(results['Primitive']["Velocity"][1:-1, -1]), label=labels[i+1], color=colors[i], ls=ls[i], lw=lw)
ax[1].set_ylabel(r'$u$')

for i,results in enumerate(resultsPickle):
    ax[2].plot(results['X Coords'][1:-1], normalize(results['Primitive']["Pressure"][1:-1, -1]), label=labels[i+1], color=colors[i], ls=ls[i], lw=lw)
ax[2].set_ylabel(r'$p$')








# Add legend only once for the figure
fig.legend(labels, loc='upper center', ncol=3, bbox_to_anchor=(0.5, +1.2))

for row in ax:
    row.set_xlabel(r'$x$')

    row.grid(alpha=.3)
plt.tight_layout()
plt.savefig(outFolder + '/flux_limiters_100.pdf', bbox_inches='tight')



plt.show()