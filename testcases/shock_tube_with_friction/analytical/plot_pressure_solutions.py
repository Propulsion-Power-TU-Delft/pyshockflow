import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
from pyshockflow import RiemannProblem
from pyshockflow.plot_styles import *

markers = ['-o', '-s', '-D', '-^', '-v']
pointInterval = 2

inputFiles = ['solutions/Test%i.pik' % i for i in range(1, 6)]
for iInput, inputFile in enumerate(inputFiles):
    with open(inputFile, 'rb') as file:
        riem = pickle.load(file)
        xCoords = riem.x        
        plt.figure()
        timeIstants = np.linspace(0, len(riem.t), 2, dtype=int, endpoint=False)
        iMark = 0
        for i in timeIstants:
            plt.plot(xCoords[::pointInterval], riem.p[::pointInterval,i], markers[iMark], mfc='none', label=r'$t=%.3f \ \rm{[s]}$' % riem.t[i])
            iMark += 1
        plt.plot(xCoords[::pointInterval], riem.p[::pointInterval,-1], markers[iMark], mfc='none', label=r'$t=%.3f \ \rm{[s]}$' % riem.t[-1])
        plt.xlabel(r'$x \ \rm{[-]}$')
        plt.ylabel(r'$p \ \rm{[-]}$')
        plt.grid(alpha=.3)
        plt.legend()
        plt.savefig('pictures/pressure_riemann_%i.pdf' %(iInput+1), bbox_inches='tight')
    
    

plt.show()
