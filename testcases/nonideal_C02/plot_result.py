import matplotlib.pyplot as plt
import numpy as np
import pickle
from pyshockflow import RiemannProblem
from pyshockflow import Driver
from pyshockflow.thesis_plots import *
import os
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset


solFile_1 = "Results/ideal_standard_NX_500/Results.pik"
solFile_3 = "Results/real_arabi_NX_500/Results.pik"

outFolder = 'Pictures'
os.makedirs(outFolder, exist_ok=True)

rhoRef = np.loadtxt("ReferenceData/rho.txt", skiprows=1)
pRef = np.loadtxt("ReferenceData/p.txt", skiprows=1)
uRef = np.loadtxt("ReferenceData/u.txt", skiprows=1)

# reference file
with open(solFile_1, 'rb') as file:
    sol1 = pickle.load(file)
with open(solFile_3, 'rb') as file:
    sol3 = pickle.load(file)
sols = [
    sol1, 
    sol3, 
    ]
labels = [
    'Standard Roe', 
    'Generalized Roe', 
    ]

set_thesis_style()
msbig = 3
ms_small = 2.0
fig, ax = create_figure(fraction=1.0, aspect_ratio=1.0, subplots=(1, 3))

ax[0].plot(rhoRef[:,0], rhoRef[:,1], 'ko', mfc='none', ms=msbig, label="Reference")
for i,res in enumerate(sols):
    ax[0].plot(res['X Coords'], res['Primitive']['Density'][:, -1], label=labels[i]) 
ax[0].set_ylabel(r'$\rho \ \rm{[kg/m^3]}$')
# Inset for Plot 1
inset1 = inset_axes(ax[0], width="35%", height="35%", loc='upper right')
inset1.plot(rhoRef[:,0], rhoRef[:,1], 'ko', mfc='none', ms=ms_small)
for i, res in enumerate(sols):
    inset1.plot(res['X Coords'], res['Primitive']['Density'][:, -1])
inset1.set_xlim(7.5, 9.0)  
inset1.set_ylim(0, 20)  
inset1.tick_params(labelsize=8)
mark_inset(ax[0], inset1, loc1=2, loc2=4, fc="none", ec="0.5")




ax[1].plot(uRef[:,0], uRef[:,1], 'ko', mfc='none', ms=msbig)
for i,res in enumerate(sols):
    ax[1].plot(res['X Coords'], res['Primitive']['Velocity'][:, -1]) 
ax[1].set_ylabel(r'$u \ \rm{[m/s]}$')
# Inset for Plot 3
# inset3 = inset_axes(ax[1], width="35%", height="35%", loc='upper left')
inset3 = ax[1].inset_axes([0.18, 0.6, 0.35, 0.35])
inset3.plot(uRef[:,0], uRef[:,1], 'ko', mfc='none', ms=ms_small)
for i, res in enumerate(sols):
    inset3.plot(res['X Coords'], res['Primitive']['Velocity'][:, -1])
inset3.set_xlim(5.7, 8.8)
inset3.set_ylim(880, 960)
mark_inset(ax[1], inset3, loc1=2, loc2=4, fc="none", ec="0.5")



ax[2].plot(pRef[:,0], pRef[:,1], 'ko', mfc='none', ms=msbig)
for i,res in enumerate(sols):
    ax[2].plot(res['X Coords'], res['Primitive']['Pressure'][:, -1]/1e6) 
ax[2].set_ylabel(r'$p \ \rm{[MPa]}$')
# Inset for Plot 2
inset2 = inset_axes(ax[2], width="35%", height="35%", loc='upper right')
inset2.plot(pRef[:,0], pRef[:,1], 'ko', mfc='none', ms=ms_small)
for i, res in enumerate(sols):
    inset2.plot(res['X Coords'], res['Primitive']['Pressure'][:, -1]/1e6)
inset2.set_xlim(6.1, 9.3)
inset2.set_ylim(0, 6)
inset2.tick_params(labelsize=8)
mark_inset(ax[2], inset2, loc1=2, loc2=4, fc="none", ec="0.5")


# for inset in [inset1, inset2, inset3]:
#     inset.set_xticks([])
#     inset.set_yticks([])


for row in ax:
    row.set_xlabel(r'$x \ \rm{[m]}$')
    row.grid(alpha=.3)
    fig.legend(loc='outside upper center', ncol=5)
# fig.legend(loc='upper center', bbox_to_anchor=(0.55, 1.1), ncol=5)
plt.savefig(outFolder + '/co2.pdf')





plt.show()