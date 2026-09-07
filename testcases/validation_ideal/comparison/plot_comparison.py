import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
from pyshockflow import RiemannProblem
from pyshockflow.thesis_plots import *

testNumber = [1,3,5]

#data
analyticalResults = ['../analytical/solutions/Test%i.pik' % i for i in testNumber]
godunovResults = ['../godunov/Results/Test%i_NX_100/Results.pik' % i for i in testNumber]
roeResults = ['../roe/Results/Test%i_NX_100/Results.pik' % i for i in testNumber]

set_thesis_style()

for iInput in range(len(analyticalResults)):
    fig, ax = create_figure(fraction=1.0, aspect_ratio=1.0, subplots=(1, 3), is_print=False)
        
    # ANALYTICAL
    with open(analyticalResults[iInput], 'rb') as file:
        res = pickle.load(file)
        ax[0].plot(res.x+0.5, res.rho[:,-1], label=r'Reference')
        ax[1].plot(res.x+0.5, res.u[:,-1])
        ax[2].plot(res.x+0.5, res.p[:,-1])
        
    # GODUNOV
    with open(godunovResults[iInput], 'rb') as file:
        res = pickle.load(file)
        ax[0].plot(
            res['X Coords'][1:-1], 
            res['Primitive']['Density'][1:-1,-1], 
            '--', 
            mfc='none', 
            label=r'Godunov')
        ax[1].plot(
            res['X Coords'][1:-1], 
            res['Primitive']['Velocity'][1:-1,-1], 
            '--', 
            mfc='none')
        ax[2].plot(
            res['X Coords'][1:-1], 
            res['Primitive']['Pressure'][1:-1,-1], 
            '--', 
            mfc='none')
    
    # ROE
    with open(roeResults[iInput], 'rb') as file:
        res = pickle.load(file)
        ax[0].plot(
            res['X Coords'][1:-1], 
            res['Primitive']['Density'][1:-1,-1],
            '-.', label=r'Roe', mfc='none')
        ax[1].plot(
            res['X Coords'][1:-1], 
            res['Primitive']['Velocity'][1:-1,-1],
            '-.')
        ax[2].plot(
            res['X Coords'][1:-1], 
            res['Primitive']['Pressure'][1:-1,-1],
            '-.')
    
    for axx in ax:
            axx.set_xlabel(r'$x$')
            axx.grid(alpha=.3)
    
    ax[0].set_ylabel(r'$\rho$')
    ax[1].set_ylabel(r'$u$')
    ax[2].set_ylabel(r'$p$')
    
    # Add a single legend at the bottom center
    fig.legend(loc='outside upper center', ncol=3)
    # Adjust layout to make room for the legend
    
    fig.savefig('Pictures/Test%i.pdf' % testNumber[iInput])


    
    

plt.show()
