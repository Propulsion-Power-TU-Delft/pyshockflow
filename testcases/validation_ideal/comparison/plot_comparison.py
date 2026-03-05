import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
from pyshockflow import RiemannProblem
from pyshockflow import *

testNumber = [1,3,5]
markerSize = 4
# ANALYTICAL SOLUTIONS
analyticalResults = ['../analytical/solutions/Test%i.pik' % i for i in testNumber]
godunovResults = ['../godunov/Results/Test%i_NX_250/Results.pik' % i for i in testNumber]
roeResults = ['../roe/Results/Test%i_NX_250/Results.pik' % i for i in testNumber]
lw=2

figSize = (10, 3.5)
for iInput in range(len(analyticalResults)):
    
    fig, ax = plt.subplots(1,3, figsize=figSize)
    
    # ANALYTICAL
    with open(analyticalResults[iInput], 'rb') as file:
        res = pickle.load(file)
        ax[0].plot(res.x+0.5, res.rho[:,-1], label=r'Reference', lw=lw)
        ax[1].plot(res.x+0.5, res.u[:,-1], lw=lw)
        ax[2].plot(res.x+0.5, res.p[:,-1], lw=lw)
        
        ax[0].set_ylabel(r'$\rho$')
        ax[1].set_ylabel(r'$u$')
        ax[2].set_ylabel(r'$p$')
        
    # GODUNOV
    with open(godunovResults[iInput], 'rb') as file:
        res = pickle.load(file)
        ax[0].plot(
            res['X Coords'][1:-1], 
            res['Primitive']['Density'][1:-1,-1], 
            '--', 
            ms=markerSize, 
            mfc='none', 
            label=r'Godunov',
            lw=lw)
        ax[1].plot(
            res['X Coords'][1:-1], 
            res['Primitive']['Velocity'][1:-1,-1], 
            '--', 
            ms=markerSize, 
            mfc='none')
        ax[2].plot(
            res['X Coords'][1:-1], 
            res['Primitive']['Pressure'][1:-1,-1], 
            '--', 
            ms=markerSize, 
            mfc='none',
            lw=lw)
    
    # ROE
    with open(roeResults[iInput], 'rb') as file:
        res = pickle.load(file)
        ax[0].plot(
            res['X Coords'][1:-1], 
            res['Primitive']['Density'][1:-1,-1],
            '-.', ms=markerSize, label=r'Roe', mfc='none', lw=lw)
        ax[1].plot(
            res['X Coords'][1:-1], 
            res['Primitive']['Velocity'][1:-1,-1],
            '-.', ms=markerSize, mfc='none', lw=lw)
        ax[2].plot(
            res['X Coords'][1:-1], 
            res['Primitive']['Pressure'][1:-1,-1],
            '-.', ms=markerSize, mfc='none', lw=lw)
        
    
    for axx in ax:
            axx.set_xlabel(r'$x$')
            axx.grid(alpha=.3)
    
    # Add a single legend at the bottom center
    # Add a single legend at the bottom center
    fig.legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, +1.11))
    plt.tight_layout()
    # Adjust layout to make room for the legend
    
    fig.savefig('Pictures/Test%i.pdf' % testNumber[iInput], bbox_inches='tight')


    
    

plt.show()
