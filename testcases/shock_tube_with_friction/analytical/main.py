import matplotlib.pyplot as plt
import numpy as np
from pyshockflow import RiemannProblem


# Page 129 of Riemann Solvers and numerical methods for fluid dynamics by Toro et al.
initialCond = {'Test1': [1.0, 0.125, 0.0, 0.0, 1.0, 0.1],
               'Test2': [1.0, 1.0, -2.0, 2.0, 0.4, 0.4],
               'Test3': [1.0, 1.0, 0.0, 0.0, 1000.0, 0.01],
               'Test4': [1.0, 1.0, 0.0, 0.0, 0.01, 100.0],
               'Test5': [5.99924, 5.99242, 19.5975, -6.19633, 460.894, 46.0950]}

nt = 5  # number of time instants on which to compute the solutions
timeVectors = {'Test1': np.linspace(0, 0.25, nt),
               'Test2': np.linspace(0, 0.15, nt),
               'Test3': np.linspace(0, 0.012, nt),
               'Test4': np.linspace(0, 0.035, nt),
               'Test5': np.linspace(0, 0.035, nt)}

nx = 250  # number of points along x
x = np.linspace(-0.5, 0.5, nx)

for key in initialCond.keys():
    initalCondition = initialCond[key]
    t = timeVectors[key]
    riem = RiemannProblem(x, t)
    riem.initializeState(initalCondition)
    riem.initializeSolutionArrays()
    riem.computeStarRegion()
    riem.solve()
    riem.showAnimation()
    riem.drawSpaceTimePlot(folder_name='pictures', file_name=key)
    riem.saveSolution(folder_name='solutions', file_name=key)
    riem.plotSolution(-1, folder_name='pictures', file_name=key)

plt.show()
