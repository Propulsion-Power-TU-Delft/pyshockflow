import CoolProp.CoolProp as CP
import numpy as np
import matplotlib.pyplot as plt
import os
import pickle
import shutil
from pathlib import Path
from pyshockflow import FluidIdeal, FluidReal

class Output():
    
    def __init__(self, filepath):
        files = [f for f in os.listdir(filepath) if os.path.isfile(os.path.join(filepath, f)) and 'pik' in f]
        files = sorted(files)
        nTimes = len(files)
        
        if nTimes == 0:
            raise ValueError('No files found in the directory')
        elif nTimes > 1: # reassemble the results in a single pickle file
            self.regroupSingleResults(filepath, files)
        else:
            self.readGlobalResult(filepath, files[0])
        
    
    def regroupSingleResults(self, filepath, files):
        nTimes = len(files)
        self.time = np.zeros(nTimes)
        self.solution = {}
        
        print("Regrouping all the results in a single file...")
        for iFile in range(len(files)):
            print(f"Reading File {iFile+1} of {len(files)}")
            with open(filepath / files[iFile], 'rb') as file:
                result = pickle.load(file)
                
                if iFile == 0:
                    self.xNodes = result['X Coords'][1:-1]
                    nNodes = len(self.xNodes)
                    self.area = result['Area Tube'][1:-1]
                    self.iterationCounter = result['Iteration Counter']
                    self.timeVec = np.zeros(nTimes)
                    self.solution['Density'] = np.zeros((nNodes, nTimes))
                    self.solution['Velocity'] = np.zeros((nNodes, nTimes))
                    self.solution['Pressure'] = np.zeros((nNodes, nTimes))
                
                self.timeVec[iFile] = result['Time']
                self.solution['Density'][:, iFile] = result['Primitive']['Density'][1:-1]
                self.solution['Velocity'][:, iFile] = result['Primitive']['Velocity'][1:-1]
                self.solution['Pressure'][:, iFile] = result['Primitive']['Pressure'][1:-1]
        
        globalOutput = {'X Coords': self.xNodes, 
                        'Area': self.area,
                        'Time': self.timeVec, 
                        'Solution': self.solution}
        
        print("Replacing all individual files with a single pickle (this could take a while) ...")
        shutil.rmtree(filepath)
        os.makedirs(filepath, exist_ok=True)
        with open(filepath / 'Results.pik', 'wb') as file:
            pickle.dump(globalOutput, file)
        print(f"Regrouped all the times in a single file: {filepath / 'Results.pik'}")
    
    
    def readGlobalResult(self, filepath, inputFile):
        with open(filepath / inputFile, 'rb') as file:
            result = pickle.load(file)
        
        self.xNodesVirtual = result['X Coords']
        self.area = result['Area']
        self.timeVec = result['Time']
        self.solution = result['Primitive']
        self.fluid = result['Fluid']
        self.config = result['Configuration']
        
        print(f"Read the file: {filepath / inputFile}")
    
    
    def showAnimation(self, jumpInstants=250):
        """
        Show animation of the results at all time instants
        """
        ni, nt = self.solution['Density'].shape
        
        def plot_limits(f, extension=0.05):
            max = f.max()
            min = f.min()
            left = min-(max-min)*extension
            right = max+(max-min)*extension
            return left, right
        
        fig, ax = plt.subplots(2, 2, figsize=(12, 8))
        density_limits = plot_limits(self.solution['Density'])
        velocity_limits = plot_limits(self.solution['Velocity'])
        pressure_limits = plot_limits(self.solution['Pressure'])
        
        # compute mach number
        if isinstance(self.fluid, FluidIdeal):
            mach = self.fluid.computeMach_u_p_rho(self.solution['Velocity'], self.solution['Pressure'], self.solution['Density'])
        elif isinstance(self.fluid, FluidReal):
            mach = np.zeros((ni, nt))
            for i in range(ni):
                mach[i, :] = self.fluid.computeMach_u_p_rho(self.solution['Velocity'][i, :], self.solution['Pressure'][i, :], self.solution['Density'][i, :])
        else:
            raise ValueError('Unknown fluid type')
        mach_limits = plot_limits(mach)
        
        interval = jumpInstants
        for it in range(0, len(self.timeVec), interval):
            for row in ax:
                for col in row:
                    col.cla()
            ax[0, 0].plot(self.xNodesVirtual[1:-1], self.solution['Density'][1:-1, it], '-C0', ms=2)
            ax[0, 0].set_ylabel(r'Density [kg/m3]')
            ax[0, 0].set_ylim(density_limits)

            ax[0, 1].plot(self.xNodesVirtual[1:-1], self.solution['Velocity'][1:-1, it], '-C1', ms=2)
            ax[0, 1].set_ylabel(r'Velocity [m/s]')
            ax[0, 1].set_ylim(velocity_limits)

            ax[1, 0].plot(self.xNodesVirtual[1:-1], self.solution['Pressure'][1:-1, it], '-C2', ms=2)
            ax[1, 0].set_ylabel(r'Pressure [Pa]')
            ax[1, 0].set_ylim(pressure_limits)

            ax[1, 1].plot(self.xNodesVirtual[1:-1], mach[1:-1, it], '-C3', ms=2)
            ax[1, 1].set_ylabel(r'Mach [-]')
            ax[1, 1].set_ylim(mach_limits)

            fig.suptitle('Time %.3e [s]' % self.timeVec[it])
            fig.tight_layout()
            for row in ax:
                for col in row:
                    col.set_xlabel('x')
                    col.grid(alpha=.3)
            plt.pause(1e-6)
        
