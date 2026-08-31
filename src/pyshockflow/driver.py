import numpy as np
import matplotlib.pyplot as plt
import os
import pickle
import csv
import sys
from pathlib import Path
from pyshockflow import RiemannProblem
from pyshockflow import AdvectionRoeBase, AdvectionRoeArabi, AdvectionRoeVinokur
from pyshockflow import FluidIdeal, FluidReal
from pyshockflow.output import Output
from pyshockflow.math_utils import *
from pyshockflow.roe_vectorized import (compute_roe_flux_ideal,
                                        compute_roe_flux_arabi,
                                        compute_roe_flux_vinokur)


class Driver:
    def __init__(self, config):
        """
        Initializes the problem with space and time arrays, along with additional fluid properties.

        Parameters
        ----------
        config: configuration file object

        Returns
        -------
        None
        """
        self.config = config
        self.topology = self.config.getTopology()
        self.isWallFrictionActive = self.config.isWallFrictionActive()
        self.frictionCoefficient = self.config.getFrictionCoefficient()
        self.isWallHeatTransferActive = self.config.isWallHeatTransferActive()
        self.wallHeatFlux = self.config.getWallHeatFlux()
        
        self.fluidName = self.config.getFluidName()
        self.fluidModel = self.config.getFluidModel()
        
        if self.fluidModel.lower()=='ideal':
            self.gmma = self.config.getFluidGamma()
            self.Rgas = self.config.getGasRConstant()
            self.fluid = FluidIdeal(self.gmma,self.Rgas)
        elif self.fluidModel.lower()=='real':
            fluidLibrary = self.config.getFluidLibrary()
            tmp = ['RefProp', 'CoolProp', 'StanMix', 'PCP-SAFT']
            if fluidLibrary not in tmp:
                raise ValueError(f"Invalid fluid library: {fluidLibrary}. Must be one of {tmp}")
            self.fluid = FluidReal(self.fluidName, fluidLibrary, False)
        
        # fluid initial states
        self.pressureLeft = self.config.getPressureLeft()
        self.pressureRight = self.config.getPressureRight()
        
        try:
            self.densityLeft = self.config.getDensityLeft()
            self.densityRight = self.config.getDensityRight()
            self.temperatureLeft = self.fluid.computeTemperature_p_rho(self.pressureLeft, self.densityLeft)
            self.temperatureRight = self.fluid.computeTemperature_p_rho(self.pressureRight, self.densityRight)
        except:
            self.temperatureLeft = self.config.getTemperatureLeft()
            self.temperatureRight = self.config.getTemperatureRight()
            self.densityLeft = self.fluid.computeDensity_p_T(self.pressureLeft, self.temperatureLeft)
            self.densityRight = self.fluid.computeDensity_p_T(self.pressureRight, self.temperatureRight)
        
        self.velocityLeft = self.config.getVelocityLeft()
        self.velocityRight = self.config.getVelocityRight()
        self.energyLeft = self.fluid.computeStaticEnergy_p_rho(self.pressureLeft, self.densityLeft)
        self.energyRight = self.fluid.computeStaticEnergy_p_rho(self.pressureRight, self.densityRight)
        
        # geometry
        self.length = self.config.getLength()
        self.nNodes = self.config.getNumberOfPoints()
        xNodes = self.generatePhysicalGeometry(self.length, self.nNodes)
        self.generateVirtualGeometry(xNodes)
        self.prepareOutputPaths()
        
        # Time related information
        self.cflMax = self.config.getCFLMax()
        self.timeMax = self.config.getTimeMax()
        
        # Boundary Conditions
        self.boundaryType = self.config.getBoundaryConditions()    
        print("Boundary Conditions Left:                    %s" %self.boundaryType[0])
        print("Boundary Conditions Right:                   %s" %self.boundaryType[1])
        print("="*80)
        
        # Print info
        print("\n" + "=" * 80)
        print(" " * 25 + "🚀  WELCOME TO PYSHOCKTUBE 🚀")
        print(" " * 18 + "Fluid Dynamics Simulation for Shock Tubes")
        print("=" * 80)
        print()  
        print("=" * 80)
        print(" "*32 + "SIMULATION DATA")
        print("Length of the domain [m]:                    %.6e" % self.length)
        print("Number of points:                            %i" % self.nNodes)
        print("Final time instant [s]:                      %.6e" % self.timeMax)
        print("Fluid name:                                  %s" % self.fluidName)
        print("Fluid treatment:                             %s" % self.fluidModel)
        if self.fluidModel.lower()=='ideal':
            print("Fluid cp/cv ratio [-]:                       %.6e" %self.gmma)
            print("Fluid gas constant [J/kgK]:                  %.6e" %self.Rgas)
        
        self.instantiatePrimitiveArrays()
        self.instantiateConservativeArrays()
        
        restartFile = self.config.getRestartFile()
        if restartFile is not None:
            self.initializeFromRestartFile(restartFile)
        else:
            self.imposeInitialConditions()
        
        self.setBoundaryConditions()
    
    
    def prepareOutputPaths(self):
        self.resultsFolder = Path(self.config.getOutputFolder())
        self.resultsFolder.mkdir(parents=True, exist_ok=True)

        self.workingDir = Path.cwd()

        # Build path safely
        resultsFilename = f"{self.config.getOutputFileName()}_NX_{self.nNodes}"
        self.resultsPath = self.resultsFolder / resultsFilename

        dum = self.resultsPath
        counter = 1

        while dum.exists():
            dum = self.resultsFolder / f"{resultsFilename}_{counter}"
            counter += 1

        self.resultsPath = dum
        self.resultsPath.mkdir(parents=True, exist_ok=True)
            
        
                

    def generatePhysicalGeometry(self, length, nodes):
        isMeshRefined = self.config.isMeshRefined()
        if isMeshRefined is False:
            xNodes = np.linspace(0, length, nodes)
        else:
            refinementCoords = self.config.getRefinementBoundaries()
            print("Mesh is refined between the two boundaries [m]: ", refinementCoords)
            
            pointsRefinement = self.config.getNumberPointsRefinement()
            totalPoints = nodes
            pointsOutside = totalPoints-pointsRefinement
            lengthUpstream = refinementCoords[0]
            lengthDownstream = length-refinementCoords[1]
            pointsUpstream = int(pointsOutside*(lengthUpstream)/(lengthUpstream+lengthDownstream))
            pointsDownstream = pointsOutside-pointsUpstream
            
            # case in which the refinement is internal
            if pointsDownstream>0 and pointsUpstream>0:
                xUpstream = np.linspace(0, refinementCoords[0], pointsUpstream+1)
                xRefinement = np.linspace(refinementCoords[0], refinementCoords[1], pointsRefinement+1)
                xDownstream = np.linspace(refinementCoords[1], length, pointsDownstream)
                if self.config.adaptMeshRefinementExtremities():
                    xUpstream = self.computeStretchedGridPoints(xUpstream, xRefinement, 'upstream')
                    xDownstream = self.computeStretchedGridPoints(xDownstream, xRefinement, 'downstream')
                xNodes = np.concatenate((xUpstream[0:-1], xRefinement[0:-1], xDownstream))
                
            elif pointsUpstream>0 and pointsDownstream==0: # the refinement finish with the end of the domain
                xUpstream = np.linspace(0, refinementCoords[0], pointsUpstream+1)
                xRefinement = np.linspace(refinementCoords[0], refinementCoords[1], pointsRefinement+1)
                if self.config.adaptMeshRefinementExtremities():
                    xUpstream = self.computeStretchedGridPoints(xUpstream, xRefinement, 'upstream')  
                xNodes = np.concatenate((xUpstream[0:-1], xRefinement))
            
            elif pointsUpstream==0 and pointsDownstream>0: # the refinement starts with the domain
                xRefinement = np.linspace(refinementCoords[0], refinementCoords[1], pointsRefinement+1)
                xDownstream = np.linspace(refinementCoords[1], length, pointsDownstream)
                if self.config.adaptMeshRefinementExtremities():
                    xDownstream = self.computeStretchedGridPoints(xDownstream, xRefinement, 'downstream')
                xNodes = np.concatenate((xRefinement[0:-1], xDownstream))
            
            else:
                raise ValueError('The refinement is ill-positioned. Please locate it internally to the domain, or at one of the extremities')
            
            self.nNodes = len(xNodes)
        
        return xNodes  
    
    
    def computeGridSpacing(self, xNodes):
        dx = np.zeros_like(xNodes)
        dx[0] = xNodes[1]-xNodes[0]
        for i in range(1,len(dx)-1):
            dx[i] = (xNodes[i+1]-xNodes[i])/2 + (xNodes[i]-xNodes[i-1])/2
        dx[-1] = xNodes[-1]-xNodes[-2]
        return dx
    
    
    def computeStretchedGridPoints(self, xCoords, xRefinement, location):      
        if location=='upstream':
            xNew = xCoords[0:-1].copy()
            dxMin = np.min(self.computeGridSpacing(xNew))
            dxRef = np.min(self.computeGridSpacing(xRefinement))
            while (dxMin-dxRef>0):
                newPoint =  xNew[-1] + (xRefinement[0]-xNew[-1])*0.5
                xNew = np.append(xNew, newPoint)
                dxMin = np.min(self.computeGridSpacing(xNew))
            xNew[-1] = xRefinement[0]
            # xNew = np.append(xNew, xRefinement[0])
        else:
            xNew = xCoords[1:].copy()
            dxMin = np.min(self.computeGridSpacing(xNew))
            dxRef = np.min(self.computeGridSpacing(xRefinement))
            while (dxMin-dxRef>0):
                newPoint =  xNew[0] - (+xNew[0] - xRefinement[-1])*0.5
                xNew = np.insert(xNew, 0, newPoint)
                dxMin = np.min(self.computeGridSpacing(xNew))
            xNew[0] = xRefinement[-1]
            # xNew = np.insert(xNew, 0, xRefinement[-1])
        return xNew
        
        
    
    def generateVirtualGeometry(self, xNodes):
        """
        Generate the virtual geometry consisting of halo nodes for boundary conditions

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        self.xNodes = xNodes
        self.nNodesHalo = self.nNodes+2
        self.xNodesVirtual = np.zeros(self.nNodesHalo)
        self.xNodesVirtual[1:-1] = self.xNodes
        self.xNodesVirtual[0] = self.xNodes[0] - (xNodes[1]-xNodes[0])
        self.xNodesVirtual[-1] = self.xNodes[-1] + (xNodes[-1]-xNodes[-2])
        self.areaReference = self.config.getAreaReference()
        self.dx = self.computeGridSpacing(self.xNodesVirtual)
        
        if self.topology.lower()=='default':
            print("The simulation proceeds with default topology: constant area")
            self.areaTube = np.zeros_like(self.xNodesVirtual)+self.areaReference
        elif self.topology.lower()=='nozzle':
            print(f"The simulation topology is: nozzle. Reading the coordinates from the nozzle file {self.config.getNozzleFilePath()}")
            self.areaTube = self.readNozzleFile(self.xNodesVirtual, self.config.getNozzleFilePath())
        else:
            raise ValueError('Unknown topology type')
        
        self.dAreaTude_dx = np.gradient(self.areaTube, self.xNodesVirtual)
        self.dV = self.areaReference * self.dx
        
    def instantiatePrimitiveArrays(self):
        """
        Instantiate the containers for the solutions. The first dimension is space, the second is time.
        """
        self.solutionNames = ['Density', 'Velocity', 'Pressure', 'Energy']
        self.solutionPrimitive = {}
        for name in self.solutionNames:
            self.solutionPrimitive[name] = np.zeros(self.nNodesHalo)


    def instantiateConservativeArrays(self):
        """
        Instantiate the containers for the solutions. The first dimension is space, the second is time.
        """
        self.solutionConsNames = ['u1', 'u2', 'u3']
        self.solutionConservative = {}
        for name in self.solutionConsNames:
            self.solutionConservative[name] = np.zeros(self.nNodesHalo)
        

    def imposeInitialConditions(self):
        """
        Initialize the conditions based on initial state, defined by right and left values
        """
        initialConditions = {'Density': np.array([self.densityLeft, self.densityRight]), 
                             'Velocity': np.array([self.velocityLeft, self.velocityRight]), 
                             'Pressure': np.array([self.pressureLeft, self.pressureRight]),
                             'Energy': np.array([self.energyLeft, self.energyRight])}
        
        for name in self.solutionNames:
            self.solutionPrimitive[name] = self.copyInitialState(initialConditions[name][0], initialConditions[name][1])
        
        print(f"Initial L/R density values [kg/m3]:          ({self.densityLeft:.6e}, {self.densityRight:.6e})")
        print(f"Initial L/R velocity values [m/s]:           ({self.velocityLeft:.6e}, {self.velocityRight:.6e})")
        print(f"Initial L/R pressure values [Pa]:            ({self.pressureLeft:.6e}, {self.pressureRight:.6e})")
        print(f"Initial L/R temperature values [K]:          ({self.temperatureLeft:.6e}, {self.temperatureRight:.6e})")
        print(f"Initial L/R energy values [J/kg]:            ({self.energyLeft:.6e}, {self.energyRight:.6e})")
        
    
    def initializeFromRestartFile(self, restartFile):
        with open(restartFile, 'rb') as file:
            restartData = pickle.load(file)
                
        for name in ['Density', 'Velocity', 'Pressure']:
            self.solutionPrimitive[name] = np.interp(self.xNodesVirtual, restartData['X Coords'], restartData['Primitive'][name][:,-1])
        
        for i in range(self.solutionPrimitive['Energy'].shape[0]):
            self.solutionPrimitive['Energy'][i] = self.fluid.computeStaticEnergy_p_rho(self.solutionPrimitive['Pressure'][i], self.solutionPrimitive['Density'][i])
            
    

    def initialConditionsArrays(self, dictIn):
        """
        Initialize the conditions based on initial state, defined by arrays
        """
        dictIn['Energy'] = dictIn['Pressure'] / (self.gmma - 1) / dictIn['Density']
        for name in self.solutionNames:
            self.solutionPrimitive[name][:, 1:-1] = dictIn[name]
    
    
    def plotGridGeometry(self, trueAspectRatio=True, pointsToJump=1, save_filename=None):
        """Plot the grid geometry. 1D tube, with thickness equal to the diameter of the tube
        """
        diameter = np.sqrt(4*self.areaTube/np.pi)
        yLower = np.zeros_like(self.xNodesVirtual)-diameter/2
        yUpper = diameter/2
        
        plt.figure()
        plt.plot(self.xNodesVirtual, yLower, 'k')
        plt.plot(self.xNodesVirtual, yUpper, 'k')
        nPointsPic = 10
        for i in range(0, len(diameter), pointsToJump):
            plt.plot(np.zeros(nPointsPic)+self.xNodesVirtual[i], np.linspace(yLower[i], yUpper[i], nPointsPic), '-k', lw=0.5)
        plt.xlabel(r'$x \ \rm{[m]}$')
        plt.ylabel(r'$r \ \rm{[m]}$')
        
        if trueAspectRatio:
            ax = plt.gca()
            ax.set_aspect('equal')
        
        if save_filename is not None:
            plt.savefig(save_filename + '.pdf', bbox_inches='tight')
        

    def copyInitialState(self, fL, fR):
        """
        Given left and right values, copy these values along the x-axis
        :param fL:
        :param fR:
        :return:
        """
        xInterface = self.config.getInterfaceLocation()
        assert(xInterface>0), f"The interface must be located within 0 and {self.length}"
        assert(xInterface<self.length), f"The interface must be located within 0 and {self.length}"
        
        f = np.zeros_like(self.xNodesVirtual)
        for i in range(len(self.xNodesVirtual)):
            if self.xNodesVirtual[i] <= xInterface:
                f[i] = fL
            else:
                f[i] = fR
        return f

    def setBoundaryConditions(self):
        """
        Set the correct boundary condition type (`reflective`, `transparent`, or `periodic`)
        """
        if self.boundaryType[0].lower()=='reflective':
            self.setReflectiveBoundaryConditions('left')
        elif self.boundaryType[0].lower()=='transparent':
            self.setTransparentBoundaryConditions('left')
        elif self.boundaryType[0].lower()=='periodic':
            self.setPeriodicBoundaryConditions('left')
        elif self.boundaryType[0].lower()=='inlet':
            self.setInletBoundaryConditions('left')
        elif self.boundaryType[0].lower()=='outlet':
            self.setOutletBoundaryConditions('left')
        else:
            raise ValueError("Unknown boundary condition type on the left")
        
        if self.boundaryType[1].lower()=='reflective':
            self.setReflectiveBoundaryConditions('right')
        elif self.boundaryType[1].lower()=='transparent':
            self.setTransparentBoundaryConditions('right')
        elif self.boundaryType[1].lower()=='periodic':
            self.setPeriodicBoundaryConditions('right')
        elif self.boundaryType[1].lower()=='outlet':
            self.setOutletBoundaryConditions('right')
        elif self.boundaryType[1].lower()=='inlet':
            self.setInletBoundaryConditions('right')
        else:
            raise ValueError("Unknown boundary condition type on the right")
        
        # update also the conservative variable arrays based on what has been done on the primitive
        self.solutionConservative['u1'], self.solutionConservative['u2'], self.solutionConservative['u3'] = (getConservativesFromPrimitives(
            self.solutionPrimitive['Density'], self.solutionPrimitive['Velocity'], self.solutionPrimitive['Pressure'], self.fluid))


    def setReflectiveBoundaryConditions(self, location):
        """
        Set reflective BC
        """
        if location=='left':
            self.solutionPrimitive['Density'][0] = self.solutionPrimitive['Density'][1]
            self.solutionPrimitive['Velocity'][0] = -self.solutionPrimitive['Velocity'][1]
            self.solutionPrimitive['Pressure'][0] = self.solutionPrimitive['Pressure'][1]
            self.solutionPrimitive['Energy'][0] = self.solutionPrimitive['Energy'][1]
        elif location=='right':
            self.solutionPrimitive['Density'][-1] = self.solutionPrimitive['Density'][-2]
            self.solutionPrimitive['Velocity'][-1] = -self.solutionPrimitive['Velocity'][-2]
            self.solutionPrimitive['Pressure'][-1] = self.solutionPrimitive['Pressure'][-2]
            self.solutionPrimitive['Energy'][-1] = self.solutionPrimitive['Energy'][-2]
        else:
            raise ValueError('Unknown location specified')
            
    

    def setTransparentBoundaryConditions(self, location):
        """
        Set transparent BC
        """
        if location=='left':
            self.solutionPrimitive['Density'][0] = self.solutionPrimitive['Density'][1]
            self.solutionPrimitive['Velocity'][0] = self.solutionPrimitive['Velocity'][1]
            self.solutionPrimitive['Pressure'][0] = self.solutionPrimitive['Pressure'][1]
            self.solutionPrimitive['Energy'][0] = self.solutionPrimitive['Energy'][1]
        elif location=='right':
            self.solutionPrimitive['Density'][-1] = self.solutionPrimitive['Density'][-2]
            self.solutionPrimitive['Velocity'][-1] = self.solutionPrimitive['Velocity'][-2]
            self.solutionPrimitive['Pressure'][-1] = self.solutionPrimitive['Pressure'][-2]
            self.solutionPrimitive['Energy'][-1] = self.solutionPrimitive['Energy'][-2]
        else:
            raise ValueError('Unknown location specified')
        
        
        
    def setPeriodicBoundaryConditions(self, location):
        """
        Set periodic BC
        """
        if location=='left':
            self.solutionPrimitive['Density'][0] = self.solutionPrimitive['Density'][-2]
            self.solutionPrimitive['Velocity'][0] = self.solutionPrimitive['Velocity'][-2]
            self.solutionPrimitive['Pressure'][0] = self.solutionPrimitive['Pressure'][-2]
            self.solutionPrimitive['Energy'][0] = self.solutionPrimitive['Energy'][-2]
        elif location=='right':
            self.solutionPrimitive['Density'][-1] = self.solutionPrimitive['Density'][1]
            self.solutionPrimitive['Velocity'][-1] = self.solutionPrimitive['Velocity'][1]
            self.solutionPrimitive['Pressure'][-1] = self.solutionPrimitive['Pressure'][1]
            self.solutionPrimitive['Energy'][-1] = self.solutionPrimitive['Energy'][1]
        else:
            raise ValueError('Unknown location specified')
    
    
    def setInletBoundaryConditions(self, location):
        """
        Set periodic BC
        """
        # handle left and right extremities with the same code
        if location=='right':
            iHalo = -1
            iInternal = -2
        elif location=='left':
            iHalo = 0
            iInternal = 1
        else:
            raise ValueError('Unknown location specified')
        
        inletConditions = self.config.getInletConditions()
        totalPressure = inletConditions[0]
        totalTemperature = inletConditions[1]
        direction = inletConditions[2]
        
        # static pressure is the only info taken from the domain
        pressure = self.solutionPrimitive['Pressure'][iInternal]
        if pressure>=totalPressure: # avoid the problems that can cause
            pressure = 0.99*totalPressure   
        density, velocity, energy = self.fluid.computeInletQuantities(pressure, totalPressure, totalTemperature, direction)
        self.solutionPrimitive['Density'][iHalo] = density
        self.solutionPrimitive['Velocity'][iHalo] = velocity
        self.solutionPrimitive['Pressure'][iHalo] = pressure
        self.solutionPrimitive['Energy'][iHalo] = energy
    
    
    def setOutletBoundaryConditions(self, location):
        """
        Set periodic BC at time
        """
        # handle left and right extremities with the same code
        if location=='right':
            iHalo = -1
            iInternal = -2
        elif location=='left':
            iHalo = 0
            iInternal = 1
        else:
            raise ValueError('Unknown location specified')
            
        machOutlet = self.fluid.computeMach_u_p_rho(self.solutionPrimitive['Velocity'][iInternal], self.solutionPrimitive['Pressure'][iInternal], self.solutionPrimitive['Density'][iInternal])        
        if machOutlet<1:
            pressure = self.config.getOutletConditions() # the pressure is the information taken from outside
            velocity = self.solutionPrimitive['Velocity'][iInternal]
            density = self.solutionPrimitive['Density'][iInternal]
            energy = self.fluid.computeStaticEnergy_p_rho(pressure, density)        
            self.solutionPrimitive['Density'][iHalo] = density
            self.solutionPrimitive['Velocity'][iHalo] = velocity
            self.solutionPrimitive['Pressure'][iHalo] = pressure
            self.solutionPrimitive['Energy'][iHalo] = energy
        else:            
            self.setTransparentBoundaryConditions(location) # the boundary is equivalent to a transparent condition
            
            


    def solve(self):
        """
        Solve the equations explicitly in time (forward Euler) using a certain advectionScheme (`Godunov`, `Roe`, `WAF`). high_order
        specifies if applying or not high order reconstruction with limiters. At the moment only type one is working -> simply
        impose high_order=True
        """
        self.entropyFixActive = self.config.isEntropyFixActive()
        self.entropyFixCoefficient = self.config.getEntropyFixCoefficient()
        advectionScheme = self.config.getNumericalScheme()
        isMusclActive = self.config.isMusclActive()
        simulationType = self.config.getSimulationType()
        outputFrequency = self.config.getOutputFrequency()
        
        print()
        print("="*80)
        print(" "*33 + "START SOLVER")
        print("Numerical flux method: %s" %(advectionScheme))
        print("MUSCL reconstruction: %s" %isMusclActive)
        print("Entropy fix active: %s" %self.entropyFixActive)
        if self.config.getFluidModel()=='real':
            print("Real Gas model, library: %s" %self.config.getFluidLibrary())
        else:
            print("Ideal Gas model")
        if self.entropyFixActive:
            print("Entropy fix coefficient: %s" %self.entropyFixCoefficient)
        print("="*80)
        print()

        # short aliases
        primitiveOld = self.solutionPrimitive.copy()
        
        # write the initial time to a solution file
        self.writeSolution(it=0, time=0)
        
        time = 0
        iTime = 1
        
        # main loop
        while time < self.timeMax:
            dt = self.computeTimeStep(primitiveOld)
            if time + dt > self.timeMax:
                dt = self.timeMax - time
            newTime = time + dt
            
            residuals = self.computeResiduals(primitiveOld, dt)
            self.updateSolution(residuals)
            
            if simulationType=='steady':
                self.printInfoResiduals(iTime, newTime, residuals)        
            else:
                print(f"Iteration: {iTime}, Progress in Time {((newTime)/self.timeMax * 100):.3f} %")
            
            self.checkSimulationStatus(dt)
            self.setBoundaryConditions()
            
            if iTime%outputFrequency==0:
                self.writeSolution(iTime, newTime)
            
            time += dt  
            iTime += 1
        
        self.writeSolution(iTime, newTime)
            
        print(" "*34 + "END SOLVER")
        print("="*80)
        print(" "*25 + "FINAL ASSEMBLY OF THE RESULTS")
        output = Output(self.resultsPath)
        print(" "*34 + "END ASSEMBLER")
        print("="*80)
    
    
    def computeResiduals(self, primitives, dt):
        availableLimiters = ['van albada', 'van leer', 'min-mod', 'superbee', 'none']
        
        limiter = self.config.getFluxLimiter()
        if limiter not in availableLimiters:
            raise ValueError(f'Limiter not recognized! Available ones are: {availableLimiters}')
        
        advectionScheme = self.config.getNumericalScheme()
        MUSCL = self.config.isMusclActive()
        
        # compute advection fluxes on every internal interface. Per unit area
        if advectionScheme.lower() == 'godunov':
            # Keep scalar path for Godunov (to be optimized separately)
            flux = np.zeros((self.nNodes+1, 3))
            for iFace in range(flux.shape[0]):
                flux[iFace, :] = self.computeFluxVector(iFace, iFace+1, primitives, dt, advectionScheme, MUSCL, limiter)
        else:
            flux = self._computeFluxVectorized(primitives, advectionScheme, MUSCL, limiter)
        
        # compute volumetric source terms, per unit volume
        source = self.computeSourceTerms(primitives)
        
        # assemble the full residual vector on every physical node (vectorized)
        dV_int = self.dV[1:-1, np.newaxis]  # (nNodes, 1) for broadcasting
        residuals = (dt / dV_int) * ((flux[:-1, :] - flux[1:, :]) * self.areaReference + source[1:-1, :] * dV_int)
        return residuals

    
    def updateSolution(self, residuals):
        self.solutionConservative['u1'][1:-1] += residuals[:,0]
        self.solutionConservative['u2'][1:-1] += residuals[:,1]
        self.solutionConservative['u3'][1:-1] += residuals[:,2]
        self.updatePrimitivesFromConservatives()
    
    def updatePrimitivesFromConservatives(self):
        self.solutionPrimitive['Density'][1:-1], self.solutionPrimitive['Velocity'][1:-1], self.solutionPrimitive['Pressure'][1:-1], self.solutionPrimitive['Energy'][1:-1] = \
                getPrimitivesFromConservatives(self.solutionConservative['u1'][1:-1], self.solutionConservative['u2'][1:-1], self.solutionConservative['u3'][1:-1], self.fluid)
        

    def computeTimeStep(self, primitive):
        velocity = primitive['Velocity'][1:-1]
        speedOfSound = np.zeros_like(velocity)
        for i in range(len(speedOfSound)):
            speedOfSound[i] = self.fluid.computeSoundSpeed_p_rho(primitive['Pressure'][i+1], primitive['Density'][i+1])
        dtMax = np.min(self.dx[1:-1] * self.cflMax / (np.abs(velocity)+speedOfSound))
        return dtMax
    
    
    def writeSolution(self, it, time):    
        
        iterationName = 'step_%06i.pik' %(it)
        fullPath = self.resultsPath / iterationName
        outputResults = {'Time': time, 
                         'Iteration Counter': it, 
                         'X Coords': self.xNodesVirtual,
                         'Area Tube': self.areaTube,
                         'Primitive': self.solutionPrimitive, 
                         'Fluid': self.fluid,
                         'Configuration': self.config}
        with open(fullPath, 'wb') as file:
            pickle.dump(outputResults, file)
    
    
    def printInfoResiduals(self, iTime, time, residuals):
        res = np.zeros(3)
        for iEq in range(3):
            res[iEq] = np.linalg.norm(residuals[:,iEq])/len(residuals[:,iEq])
            if res[iEq]!=0:
                res[iEq] = np.log10(res[iEq])
        timeProgress = time/self.timeMax * 100
        print('Iteration %i    Progress in Time %.3f%%    Residuals: %.6f, %.6f, %.6f' %(iTime, timeProgress, res[0], res[1], res[2]))
    
    def computeSourceTerms(self, primitive):
        """compute global volumetric source terms (per unit length)
        """
        source = np.zeros((self.nNodesHalo,3))
        
        # compute the source terms if needed

        if self.topology.lower()=='nozzle':
            source += self.computeGeometricalSourceTerms(primitive)
        
        if self.isWallFrictionActive:
            source += self.computeFrictionSourceTerms(primitive)
        
        if self.isWallHeatTransferActive:
            source += self.computeHeatTransferSourceTerms(primitive)
            
        return source
    
    def computeGeometricalSourceTerms(self, primitive):
        """compute source terms related to area variations along the tube due to a nozzle. Source terms taken from 'On the numerical simulation
        of non-classical quasi-1D steady nozzle flows: Capturing sonic shocks' by Vimercati and Guardone.
        """
        totalEnergy = primitive['Energy'][:] + 0.5*primitive['Velocity']**2
        source = np.zeros((self.nNodesHalo,3))
        source[:,0] = - primitive['Density'] * primitive['Velocity']*self.dAreaTude_dx/self.areaTube
        source[:,1] = - (primitive['Density'] * primitive['Velocity']**2)*self.dAreaTude_dx/self.areaTube
        source[:,2] = - primitive['Velocity'] *(primitive['Density']*totalEnergy + primitive['Pressure'])*self.dAreaTude_dx/self.areaTube
        return source
    
    def computeFrictionSourceTerms(self, primitive):
        """compute source terms related to friction along the walls
        """
        frictionCoefficient = self.frictionCoefficient 
        hydraulicDiameter = 4*self.areaTube/(np.pi*np.sqrt(self.areaTube/np.pi)*2)
        Sm = -2.0 * frictionCoefficient * primitive['Density'] * primitive['Velocity'] * np.abs(primitive['Velocity']) / hydraulicDiameter
        
        source = np.zeros((self.nNodesHalo,3))
        source[:,1] = Sm
        return source
    
    def computeHeatTransferSourceTerms(self, primitive):
        """compute source term of energy eq. related to heat transfer along the walls.
        Must be per unit length
        """
        radius = np.sqrt(self.areaTube/np.pi)
        transferArea = np.pi*radius*2 * self.dx
        heatTransfer = transferArea * self.wallHeatFlux
        volumetricHeatTransfer = heatTransfer/self.dV
        
        source = np.zeros((self.nNodesHalo,3))
        source[:,2] = volumetricHeatTransfer
        return source
    
    
    def checkSimulationStatus(self, dt):
        """
        Check if nans or infs are detected and in that case stop the simulation and provide explanation
        """
        if np.any(np.isnan(self.solutionPrimitive['Density'])) or np.any(np.isinf(self.solutionPrimitive['Density'])) or \
            np.any(np.isnan(self.solutionPrimitive['Pressure'])) or np.any(np.isinf(self.solutionPrimitive['Pressure'])):
            print()
            print()
            print("######################  SIMULATION DIVERGED ############################")
            print('NaNs detected in density. Simulation stopped.')
            cfl = self.computeMaxCFL(dt) # use the previous time step to compute where the solution had CFL related problems
            print("Maximum CFL number found: %.3f" %(np.max(cfl)))
            print("At location x: %.3f [m]" %(self.xNodesVirtual[np.argmax(cfl)]))
            print("Visualize the plot to understand critical locations, and decrease CFL_MAX input setting.")
            print("###############################  EXIT ##################################")
            print()
            
            plt.figure()
            plt.plot(self.xNodes, cfl)
            plt.xlabel('x [m]')
            plt.ylabel('CFL [-]')
            plt.grid(alpha=.3)
            plt.show()
            sys.exit()
    
    
    def computeMaxCFL(self, dt):
        pressure = self.solutionPrimitive['Pressure'][1:-1]
        density = self.solutionPrimitive['Density'][1:-1]
        velocity = self.solutionPrimitive['Velocity'][1:-1]
        dx = self.dx[1:-1]
        soundSpeed = np.zeros_like(pressure)
        for i in range(len(soundSpeed)):
            soundSpeed = self.fluid.computeSoundSpeed_p_rho(pressure[i], density[i])
        cfl = (np.abs(velocity)+soundSpeed)*dt/dx
        return cfl
        

    def computeFluxVector(self, il, ir, primitive, dt, advectionScheme, MUSCL, limiter):
        """
        compute the flux vector at the interface between grid points `il` and `ir`, using a certain `advectionScheme`.
        """
        
        # flow reconstruction if high_order=True
        if (MUSCL and il>2 and ir<self.nNodesHalo-3):
            rhoL, uL, pL, rhoR, uR, pR = self.computeMusclReconstruction(il, ir, limiter)
        else:
            rhoL = primitive['Density'][il]
            rhoR = primitive['Density'][ir]
            uL = primitive['Velocity'][il]
            uR = primitive['Velocity'][ir]
            pL = primitive['Pressure'][il]
            pR = primitive['Pressure'][ir]            
        
        # flux calculation
        if advectionScheme.lower()=='godunov':
            if self.fluidModel!='ideal':
                raise ValueError('Godunov scheme is available only for ideal gas model')
            else:
                # Godunov flux calculation
                nx, nt = 51, 51 
                x = np.linspace(-self.dx[il]/2, self.dx[ir]/2, nx)
                t = np.linspace(0, dt, nt)
                riem = RiemannProblem(x, t)
                riem.initializeState([rhoL, rhoR, uL, uR, pL, pR])
                riem.initializeSolutionArrays()
                riem.computeStarRegion()
                riem.solve(space_domain='interface', time_domain='global') # compute Riemann solution only at x=0, but on all time instants
                rho, u, p = riem.getSolutionInTime()
                u1, u2, u3 = getConservativesFromPrimitives(rho, u, p, self.fluid)
                u1AVG, u2AVG, u3AVG = np.sum(u1)/len(u1), np.sum(u2)/len(u2), np.sum(u3)/len(u3)
                flux = computeAdvectionFluxFromConservatives(u1AVG, u2AVG, u3AVG, self.fluid) 
        elif advectionScheme.lower()=='roe':
            if self.fluidModel=='real':
                raise ValueError('Basic Roe scheme is not available for real gas model. Select Roe_Arabi or Roe_Vinokur, depending on the Roe Avg procedure that you want.')
            else:
                roe = AdvectionRoeBase(rhoL, rhoR, uL, uR, pL, pR, self.fluid)
                flux = roe.computeFlux(entropyFixActive=self.entropyFixActive, fixCoefficient=self.entropyFixCoefficient)
        elif advectionScheme.lower()=='roe_arabi':
            if self.fluidModel=='ideal':
                raise ValueError('Roe_Arabi scheme is not available for ideal gas model. Select Standard Roe scheme.')
            else:
                roe = AdvectionRoeArabi(rhoL, rhoR, uL, uR, pL, pR, self.fluid)
                flux = roe.computeFlux(entropyFixActive=self.entropyFixActive, fixCoefficient=self.entropyFixCoefficient)
        elif advectionScheme.lower()=='roe_vinokur':
                roe = AdvectionRoeVinokur(rhoL, rhoR, uL, uR, pL, pR, self.fluid)
                roe.computeAveragedVariables()
                flux = roe.computeFlux(entropyFixActive=self.entropyFixActive, fixCoefficient=self.entropyFixCoefficient)
        else:
            raise ValueError('Unknown flux method')
        
        return flux
    
    def computeMusclReconstruction(self, il, ir, limiter):
        """
        MUSCL reconstruction coupled with a certain limiter
        """
        # states left, left minus 1, right, right plus one
        U_lm = np.array([self.solutionPrimitive['Density'][il-1], self.solutionPrimitive['Velocity'][il-1], self.solutionPrimitive['Pressure'][il-1]])
        U_l = np.array([self.solutionPrimitive['Density'][il], self.solutionPrimitive['Velocity'][il], self.solutionPrimitive['Pressure'][il]])
        U_r = np.array([self.solutionPrimitive['Density'][ir], self.solutionPrimitive['Velocity'][ir], self.solutionPrimitive['Pressure'][ir]])
        U_rp = np.array([self.solutionPrimitive['Density'][ir+1], self.solutionPrimitive['Velocity'][ir+1], self.solutionPrimitive['Pressure'][ir+1]])
        
        dx_left_leftm = self.xNodes[il]-self.xNodes[il-1] # dx is always the same for now
        dx_right_left = self.xNodes[ir]-self.xNodes[il]
        dx_rightp_right = self.xNodes[ir+1]-self.xNodes[ir]
        
        # compute the smoothness indicators
        smoothnessLeft = self.computeSmoothnessIndicators(U_lm, U_l, U_r, dx_left_leftm, dx_right_left)
        smoothnessRight = self.computeSmoothnessIndicators(U_l, U_r, U_rp, dx_right_left, dx_rightp_right)
        
        # compute left and right flux limiters
        psi_left = self.computeFluxLimiter(smoothnessLeft, limiter)
        psi_right = self.computeFluxLimiter(smoothnessRight, limiter)
        
        # reconstruct left and right states
        U_l_rec = U_l+0.5*psi_left*(U_r-U_l)
        U_r_rec = U_r-0.5*psi_right*(U_rp-U_r)

        return U_l_rec[0], U_l_rec[1], U_l_rec[2], U_r_rec[0], U_r_rec[1], U_r_rec[2]


    def computeSmoothnessIndicators(self, U_left, U_central, U_right, dx_left, dx_right):
        """
        compute the array of smoothness indicators for the following flux limiter evaluation
        """
        rVector = ((U_central-U_left)/dx_left) / ((U_right-U_central)/dx_right + 1e-6)
        return rVector
    
    
    def saveSolution(self):
        """
        Save the full object as a pickle for later use
        """
        folder_name = self.config.getOutputFolder()
        os.makedirs(folder_name, exist_ok=True)
        file_name = self.config.getOutputFileName()
        full_path = folder_name+'/'+file_name+'_NX_%i_TMAX_%.6f.pik' %(self.nNodes, self.timeMax)
        with open(full_path, 'wb') as file:
            pickle.dump(self, file)
        print('Pickle object with full solution saved to ' + full_path + ' !')


    def saveNodeSolutionToCSV(self, iNode, timeInstants, folder_name, file_name):
        """
        Save the array of fluid flow quantities (P,T,s,Mach,Gamma) from the solution to a CSV file.
        """
        filePath = folder_name + '/' + file_name + '.dat'

        pressure = self.solutionPrimitive['Pressure'][iNode, :]  # Extract the pressure data (1D array)
        density = self.solutionPrimitive['Density'][iNode, :]  # Extract the density data (1D array)
        temperature = self.fluid.computeTemperature_p_rho(pressure, density)
        entropy = self.fluid.computeEntropy_p_rho(pressure, density)
        fundDerGasDynamics = self.fluid.computeFunDerGamma_p_rho(pressure, density)
        compressibilityFactor = self.fluid.computeComprFactorZ_p_rho(pressure, density)

        with open(filePath, 'w', newline='', encoding='utf-8') as file:
            writer = csv.writer(file)
            for value in range(len(timeInstants)):
                writer.writerow([timeInstants[value], pressure[value], temperature[value], density[value],
                                 entropy[value], fundDerGasDynamics[value], compressibilityFactor[value]])

        print(f"Fluid flow quantities (P,T,D,s,Gamma,Z) saved to {filePath}!")


    def computeFluxLimiter(self, r_vec, limiter):
        """
        compute the flux limiter functions.
        """
        psi = np.zeros(3)
        for i in range(len(r_vec)):
            r = r_vec[i]

            if limiter.lower() == 'van albada':
                psi[i] = (r**2+r)/(1+r**2)

            elif limiter.lower() == 'van leer':
                psi[i] = (r+np.abs(r))/(1+np.abs(r))

            elif limiter.lower() == 'min-mod':
                psi[i] = np.maximum(0, np.minimum(1, r))

            elif limiter.lower() == 'superbee':
                psi[i] = np.max(np.array([0, np.minimum(2 * r, 1), np.minimum(r, 2)]))

            elif limiter.lower() == 'none':
                psi[i] = 1 
            else:
                raise ValueError('Limiter not recognized!')
            
        return psi
    

    # ------------------------------------------------------------------
    # Vectorized flux computation (Phase 1 optimization)
    # ------------------------------------------------------------------

    def _computeFluxVectorized(self, primitives, advectionScheme, MUSCL, limiter):
        """Compute all interface fluxes at once using vectorized Roe solvers.

        Replaces the per-interface Python ``for`` loop and per-call class
        instantiation with a single call to a vectorized NumPy function.
        """
        rho = primitives['Density']    # shape (nNodesHalo,)
        u   = primitives['Velocity']
        p   = primitives['Pressure']

        nFaces = self.nNodes + 1

        # --- Build left / right states for every interface (first-order) ---
        rhoL = rho[:-1].copy()
        rhoR = rho[1:].copy()
        uL   = u[:-1].copy()
        uR   = u[1:].copy()
        pL   = p[:-1].copy()
        pR   = p[1:].copy()

        # --- Apply MUSCL reconstruction to interior faces ----------------
        if MUSCL:
            self._applyMusclReconstructionVectorized(
                primitives, rhoL, rhoR, uL, uR, pL, pR, limiter)

        # --- Dispatch to the appropriate vectorized solver ----------------
        scheme = advectionScheme.lower()

        if scheme == 'roe':
            if self.fluidModel == 'real':
                raise ValueError(
                    'Basic Roe scheme is not available for real gas model. '
                    'Select Roe_Arabi or Roe_Vinokur, depending on the Roe '
                    'Avg procedure that you want.')
            return compute_roe_flux_ideal(
                rhoL, rhoR, uL, uR, pL, pR,
                self.fluid.gmma,
                self.entropyFixActive, self.entropyFixCoefficient)

        elif scheme == 'roe_arabi':
            if self.fluidModel == 'ideal':
                raise ValueError(
                    'Roe_Arabi scheme is not available for ideal gas model. '
                    'Select Standard Roe scheme.')
            # Per-interface EOS calls (scalar; real-gas vectorization in Phase 4)
            eL_arr = np.empty(nFaces)
            eR_arr = np.empty(nFaces)
            aL_arr = np.empty(nFaces)
            aR_arr = np.empty(nFaces)
            for i in range(nFaces):
                eL_arr[i] = self.fluid.computeStaticEnergy_p_rho(pL[i], rhoL[i])
                eR_arr[i] = self.fluid.computeStaticEnergy_p_rho(pR[i], rhoR[i])
                aL_arr[i] = self.fluid.computeSoundSpeed_p_rho(pL[i], rhoL[i])
                aR_arr[i] = self.fluid.computeSoundSpeed_p_rho(pR[i], rhoR[i])
            return compute_roe_flux_arabi(
                rhoL, rhoR, uL, uR, pL, pR,
                eL_arr, eR_arr, aL_arr, aR_arr,
                self.entropyFixActive, self.entropyFixCoefficient)

        elif scheme == 'roe_vinokur':
            if self.fluidModel == 'ideal':
                # Ideal-gas: chi=0, kappa=gamma-1 everywhere (fully vectorized)
                gm1 = self.fluid.gmma - 1.0
                eL_arr = pL / (gm1 * rhoL)
                eR_arr = pR / (gm1 * rhoR)
                chiL_arr = np.zeros(nFaces)
                kappaL_arr = np.full(nFaces, gm1)
                chiR_arr = np.zeros(nFaces)
                kappaR_arr = np.full(nFaces, gm1)
                chiM_arr = np.zeros(nFaces)
                kappaM_arr = np.full(nFaces, gm1)
            else:
                # Real-gas: scalar EOS calls (to be optimized in Phase 4)
                eL_arr = np.empty(nFaces)
                eR_arr = np.empty(nFaces)
                chiL_arr = np.empty(nFaces)
                kappaL_arr = np.empty(nFaces)
                chiR_arr = np.empty(nFaces)
                kappaR_arr = np.empty(nFaces)
                chiM_arr = np.empty(nFaces)
                kappaM_arr = np.empty(nFaces)
                p_mean = 0.5 * (pL + pR)
                rho_mean = 0.5 * (rhoL + rhoR)
                for i in range(nFaces):
                    eL_arr[i] = self.fluid.computeStaticEnergy_p_rho(pL[i], rhoL[i])
                    eR_arr[i] = self.fluid.computeStaticEnergy_p_rho(pR[i], rhoR[i])
                    chiL_arr[i], kappaL_arr[i] = self.fluid.computeChiKappa_VinokurScheme_p_rho(pL[i], rhoL[i])
                    chiR_arr[i], kappaR_arr[i] = self.fluid.computeChiKappa_VinokurScheme_p_rho(pR[i], rhoR[i])
                    chiM_arr[i], kappaM_arr[i] = self.fluid.computeChiKappa_VinokurScheme_p_rho(p_mean[i], rho_mean[i])

            return compute_roe_flux_vinokur(
                rhoL, rhoR, uL, uR, pL, pR,
                eL_arr, eR_arr,
                chiL_arr, chiR_arr, chiM_arr,
                kappaL_arr, kappaR_arr, kappaM_arr,
                self.entropyFixActive, self.entropyFixCoefficient)
        else:
            raise ValueError(f'Unknown flux method: {advectionScheme}')


    def _applyMusclReconstructionVectorized(self, primitives, rhoL, rhoR,
                                            uL, uR, pL, pR, limiter):
        """Apply MUSCL reconstruction to interior faces, modifying L/R arrays in-place.

        The original per-face loop with temporary 3-element allocations is
        replaced by slicing the full primitive arrays into stencil arrays and
        computing all limiters at once.
        """
        rho = primitives['Density']
        u   = primitives['Velocity']
        p   = primitives['Pressure']
        # xNodes uses the physical grid (same indexing as original code)
        x = self.xNodes

        # MUSCL condition: il > 2 and ir < nNodesHalo - 3
        # => face index iFace in [3, nNodes - 3]
        ms = 3                       # first MUSCL face
        me = self.nNodes - 3         # last MUSCL face (inclusive)

        if ms > me:
            return  # grid too small for MUSCL stencil

        nMuscl = me - ms + 1

        # Stencil slices into the halo-inclusive arrays (size nNodesHalo)
        # face iFace => il=iFace, ir=iFace+1
        # stencil: il-1, il, ir, ir+1 = iFace-1 .. iFace+2
        s0 = ms - 1                  # start of il-1 range
        U_lm = np.column_stack([rho[s0       : s0 + nMuscl],
                                u[s0         : s0 + nMuscl],
                                p[s0         : s0 + nMuscl]])
        U_l  = np.column_stack([rho[s0 + 1   : s0 + 1 + nMuscl],
                                u[s0 + 1     : s0 + 1 + nMuscl],
                                p[s0 + 1     : s0 + 1 + nMuscl]])
        U_r  = np.column_stack([rho[s0 + 2   : s0 + 2 + nMuscl],
                                u[s0 + 2     : s0 + 2 + nMuscl],
                                p[s0 + 2     : s0 + 2 + nMuscl]])
        U_rp = np.column_stack([rho[s0 + 3   : s0 + 3 + nMuscl],
                                u[s0 + 3     : s0 + 3 + nMuscl],
                                p[s0 + 3     : s0 + 3 + nMuscl]])

        # Grid spacing (uses self.xNodes — matches the original code)
        dx_lm_l = x[ms     : me + 1] - x[ms - 1 : me]        # x[il]-x[il-1]
        dx_l_r  = x[ms + 1 : me + 2] - x[ms     : me + 1]    # x[ir]-x[il]
        dx_r_rp = x[ms + 2 : me + 3] - x[ms + 1 : me + 2]    # x[ir+1]-x[ir]

        # Smoothness indicators  (nMuscl, 3)
        r_left  = ((U_l - U_lm) / dx_lm_l[:, np.newaxis]) / \
                  ((U_r - U_l)  / dx_l_r[:, np.newaxis] + 1e-6)
        r_right = ((U_r - U_l)  / dx_l_r[:, np.newaxis]) / \
                  ((U_rp - U_r) / dx_r_rp[:, np.newaxis] + 1e-6)

        # Flux limiters (nMuscl, 3)
        psi_left  = self._computeFluxLimiterVec(r_left, limiter)
        psi_right = self._computeFluxLimiterVec(r_right, limiter)

        # Reconstructed states
        U_l_rec = U_l + 0.5 * psi_left  * (U_r  - U_l)
        U_r_rec = U_r - 0.5 * psi_right * (U_rp - U_r)

        # Overwrite the L/R arrays for MUSCL faces
        rhoL[ms : me + 1] = U_l_rec[:, 0]
        uL  [ms : me + 1] = U_l_rec[:, 1]
        pL  [ms : me + 1] = U_l_rec[:, 2]
        rhoR[ms : me + 1] = U_r_rec[:, 0]
        uR  [ms : me + 1] = U_r_rec[:, 1]
        pR  [ms : me + 1] = U_r_rec[:, 2]


    @staticmethod
    def _computeFluxLimiterVec(r, limiter):
        """Vectorized flux limiter for arrays of shape (N, 3).

        The limiter type is resolved once (no per-element string comparison).
        """
        lim = limiter.lower()
        if lim == 'van albada':
            return (r ** 2 + r) / (1.0 + r ** 2)
        elif lim == 'van leer':
            abs_r = np.abs(r)
            return (r + abs_r) / (1.0 + abs_r)
        elif lim == 'min-mod':
            return np.maximum(0.0, np.minimum(1.0, r))
        elif lim == 'superbee':
            return np.maximum(0.0, np.maximum(np.minimum(2.0 * r, 1.0),
                                              np.minimum(r, 2.0)))
        elif lim == 'none':
            return np.ones_like(r)
        else:
            raise ValueError(f'Limiter not recognized: {limiter}')

    
    def readNozzleFile(self, xTube, filepath):
        nozzleData = np.loadtxt(filepath, skiprows=1, delimiter=',', dtype=float)
        nozzleX = nozzleData[:,0]
        nozzleArea = nozzleData[:,1]
        
        # Linear interpolation with external filling set to area Reference (=Tube area)
        areaReference = self.config.getAreaReference()
        interpolatedNozzleArea = np.interp(xTube, nozzleX, nozzleArea, left=areaReference, right=areaReference)
    
        print(f"The reference tube area is: {areaReference:.6f} [m2].")
        print(f"The nozzle throat area is {interpolatedNozzleArea.min():.6f} [m2].")
        print(f"The nozzle maximum area is {interpolatedNozzleArea.max():.6f} [m2].")
        print(f"The area ratio between nozzle throat and exit section is {interpolatedNozzleArea.min()/interpolatedNozzleArea[-1]:.6f}.")
        print(f"The area ratio between nozzle throat and tube is {interpolatedNozzleArea.min()/areaReference:.6f}.")
        print(f"If this is not correct, modify the REFERENCE_AREA setting in the geometry section of the input file to the correct value for the tube area, or modify the nozzle csv file to be consistent with the tube area.")
        
        return interpolatedNozzleArea
    























