import matplotlib.pyplot as plt
import numpy as np
import pickle
from pyshockflow import RiemannProblem
from pyshockflow import Driver
import os
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset



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
ax[1].plot(data.x+0.5, data.u[:,-1], '-k', label='Reference', lw=1.7)
ax[2].plot(data.x+0.5, data.p[:,-1], '-k', label='Reference', lw=1.7)



colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5']
ls = ['-', '--', ':', '-.', '--', '-.']

for i,results in enumerate(resultsPickle):
    ax[0].plot(
        results['X Coords'][1:-1], 
        normalize(results['Primitive']["Density"][1:-1, -1]), 
        label=labels[i+1], color=colors[i], ls=ls[i], lw=lw)
    
    ax[1].plot(
        results['X Coords'][1:-1], 
        normalize(results['Primitive']["Velocity"][1:-1, -1]), 
        label=labels[i+1], color=colors[i], ls=ls[i], lw=lw)
    
    ax[2].plot(results['X Coords'][1:-1], 
               normalize(results['Primitive']["Pressure"][1:-1, -1]), 
               label=labels[i+1], color=colors[i], ls=ls[i], lw=lw)
    
ax[0].set_ylabel(r'$\rho$')
ax[1].set_ylabel(r'$u$')
ax[2].set_ylabel(r'$p$')
    

# Insets for Plot 0
inset0 = inset_axes(ax[0], width="40%", height="40%", loc='upper right')
inset0.plot(data.x+0.5, data.rho[:,-1], '-k', lw=lw)
for i, results in enumerate(resultsPickle):
    inset0.plot(results['X Coords'][1:-1], results['Primitive']['Density'][1:-1, -1], lw=lw, color=colors[i], ls=ls[i])
inset0.set_xlim(0.6, 0.75) 
inset0.set_ylim(0.33, 0.45)
inset0.tick_params(labelsize=8)
mark_inset(ax[0], inset0, loc1=2, loc2=4, fc="none", ec="0.5")

# Insets for Plot 1
inset1 = inset_axes(ax[1], width="40%", height="40%",
                    bbox_to_anchor=(0.4, 0.25, 1.0, 1.0),  # (x, y, width, height) in axes coords [0-1]
                    bbox_transform=ax[1].transAxes,
                    loc='lower left')
inset1.plot(data.x+0.5, data.u[:,-1], '-k', lw=lw)
for i, results in enumerate(resultsPickle):
    inset1.plot(results['X Coords'][1:-1], results['Primitive']['Velocity'][1:-1, -1], lw=lw, color=colors[i], ls=ls[i])
inset1.set_xlim(0.806, 0.96) 
inset1.set_ylim(0.806, 0.952)
inset1.tick_params(labelsize=8)
mark_inset(ax[1], inset1, loc1=2, loc2=4, fc="none", ec="0.5")

# # Insets for Plot 2
inset2 = inset_axes(ax[2], width="40%", height="40%", loc='upper right')
inset2.plot(data.x+0.5, data.p[:,-1], '-k', lw=lw)
for i, results in enumerate(resultsPickle):
    inset2.plot(results['X Coords'][1:-1], results['Primitive']['Pressure'][1:-1, -1], lw=lw, color=colors[i], ls=ls[i])
inset2.set_xlim(0.85, 0.986) 
inset2.set_ylim(0.08, 0.32)
inset2.tick_params(labelsize=8)
mark_inset(ax[2], inset2, loc1=2, loc2=4, fc="none", ec="0.5")






# Add legend only once for the figure
fig.legend(labels, loc='upper center', ncol=3, bbox_to_anchor=(0.5, +1.2))

for row in ax:
    row.set_xlabel(r'$x$')

    row.grid(alpha=.3)

for inset in [inset0, inset1, inset2]:
    inset.set_xticks([])
    inset.set_yticks([])
    
plt.tight_layout()
plt.savefig(outFolder + '/flux_limiters_100.pdf', bbox_inches='tight')



plt.show()