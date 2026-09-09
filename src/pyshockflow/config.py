import configparser
import os

class Config:
    def __init__(self, config_file='input.ini'):
        if not os.path.exists(config_file):
            raise FileNotFoundError(f"Config file '{config_file}' not found.")
        
        self.config_parser = configparser.ConfigParser()
        self.config_parser.read(config_file)
    
    def getNumberOfPoints(self):
        return int(self.config_parser.get('SIMULATION', 'NUMBER_POINTS')) 
    
    def getLength(self):
        return float(self.config_parser.get('GEOMETRY', 'LENGTH')) 
    
    def getPressureLeft(self):
        return float(self.config_parser.get('SIMULATION', 'PRESSURE_LEFT')) 
    
    def getPressureRight(self):
        return float(self.config_parser.get('SIMULATION', 'PRESSURE_RIGHT')) 
    
    def getDensityLeft(self):
        return float(self.config_parser.get('SIMULATION', 'DENSITY_LEFT')) 
    
    def getDensityRight(self):
        return float(self.config_parser.get('SIMULATION', 'DENSITY_RIGHT')) 
    
    def getTemperatureLeft(self):
        return float(self.config_parser.get('SIMULATION', 'TEMPERATURE_LEFT')) 
    
    def getTemperatureRight(self):
        return float(self.config_parser.get('SIMULATION', 'TEMPERATURE_RIGHT')) 
    
    def getVelocityLeft(self):
        return float(self.config_parser.get('SIMULATION', 'VELOCITY_LEFT')) 
    
    def getVelocityRight(self):
        return float(self.config_parser.get('SIMULATION', 'VELOCITY_RIGHT')) 
    
    def getCFLMax(self):
        return float(self.config_parser.get('SIMULATION', 'CFL_MAX')) 
    
    def getTimeMax(self):
        return float(self.config_parser.get('SIMULATION', 'TIME_MAX')) 
    
    def getTimeStepMethod(self):
        try:
            return str(self.config_parser.get('SIMULATION', 'TIME_STEP_METHOD')).lower()
        except:
            return 'constant'
    
    def getFluidName(self):
        return str(self.config_parser.get('FLUID', 'FLUID_NAME'))
    
    def getFluidModel(self):
        return str(self.config_parser.get('FLUID', 'FLUID_MODEL')).lower() 
    
    def getFluidGamma(self):
        return float(self.config_parser.get('FLUID', 'FLUID_GAMMA'))
    
    def getGasRConstant(self):
        return float(self.config_parser.get('FLUID', 'GAS_R_CONSTANT'))
    
    def getInterfaceLocation(self):
        return float(self.config_parser.get('GEOMETRY', 'INTERFACE_LOCATION'))
    
    def getBoundaryConditions(self):
        left = str(self.config_parser.get('SIMULATION', 'BOUNDARY_CONDITION_LEFT')).lower() 
        right = str(self.config_parser.get('SIMULATION', 'BOUNDARY_CONDITION_RIGHT')).lower() 
        return (left, right)

    def getInletConditions(self):
        res = str(self.config_parser.get('SIMULATION', 'INLET_CONDITIONS')).lower() 
        res = [float(value.strip()) for value in res.split(',')]
        return res
    
    def getOutletConditions(self):
        res = float(self.config_parser.get('SIMULATION', 'OUTLET_CONDITIONS'))
        return res
    
    def getNumericalScheme(self):
        return str(self.config_parser.get('SIMULATION', 'NUMERICAL_SCHEME')).lower() 
    
    def getFluxLimiter(self):
        try:
            return str(self.config_parser.get('SIMULATION', 'FLUX_LIMITER')).lower() 
        except:
            return 'van albada'
    
    def getOutputFolder(self):
        return str(self.config_parser.get('OUTPUT', 'FOLDER_NAME')) 
    
    def getOutputFileName(self):
        return str(self.config_parser.get('OUTPUT', 'FILE_NAME')) 
    
    def showAnimation(self):
        res = str(self.config_parser.get('OUTPUT', 'SHOW_ANIMATION')).lower() 
        if res=='yes' or res=='true':
            return True
        else:
            return False
    
    
    def getRestartFile(self):
        try:
            return str(self.config_parser.get('SIMULATION', 'RESTART_FILE'))
        except:
            return None
    
    def getSimulationType(self):
        try:
            return str(self.config_parser.get('SIMULATION', 'SIMULATION_TYPE')).lower()
        except:
            return "unsteady"
    
    def isMusclActive(self):
        try:
            res = str(self.config_parser.get('SIMULATION', 'MUSCL_RECONSTRUCTION')).lower() 
            if res=='yes' or res=='true':
                return True
            else:
                return False
        except:
            return False # false by default
    
    
    def isMeshRefined(self):
        try:
            res = str(self.config_parser.get('SIMULATION', 'MESH_REFINEMENT')).lower() 
            if res=='yes' or res=='true':
                return True
            else:
                return False
        except:
            return False # false by default
    
    def isEntropyFixActive(self):
        try:
            res = str(self.config_parser.get('SIMULATION', 'ENTROPY_FIX_ACTIVE')).lower() 
            if res=='yes' or res=='true':
                return True
            else:
                return False
        except:
            return True # true by default
    
    def getEntropyFixCoefficient(self):
        try:
            res = float(self.config_parser.get('SIMULATION', 'ENTROPY_FIX_COEFFICIENT'))
            return res
        except:
            return 0.2 # 0.2 by default
    
    def getFluidLibrary(self):
        return str(self.config_parser.get('FLUID', 'FLUID_LIBRARY'))
    
    def adaptMeshRefinementExtremities(self):
        try:
            res = str(self.config_parser.get('SIMULATION', 'ADAPT_MESH_REFINEMENT')).lower() 
            if res=='yes' or res=='true':
                return True
            else:
                return False
        except:
            return False #  by default
    
    def getRefinementBoundaries(self):
        start = float(self.config_parser.get('SIMULATION', 'X_START_REFINEMENT')) 
        end = float(self.config_parser.get('SIMULATION', 'X_END_REFINEMENT')) 
        return (start, end)
    
    def getNumberPointsRefinement(self):
        return int(self.config_parser.get('SIMULATION', 'NUMBER_POINTS_REFINEMENT')) 
        
    
    def getTopology(self):
        try:
            return str(self.config_parser.get('GEOMETRY', 'TOPOLOGY')).lower()
        except:
            return 'default' 
    
    def getFrictionCoefficient(self):
        try:
            return float(self.config_parser.get('SIMULATION', 'FRICTION_COEFFICIENT'))
        except:
            return 0.003 # default value
        
    def isWallFrictionActive(self):
        try:
            res = str(self.config_parser.get('SIMULATION', 'WALL_FRICTION_ACTIVE')).lower()
            if res=='yes' or res=='true':
                return True
            else:
                return False
        except:
            return False 

    def getWallFrictionModel(self):
        try:
            return str(self.config_parser.get('SIMULATION', 'WALL_FRICTION_MODEL')).lower()
        except:
            return 'constant'

    def getFrictionTransitionReynolds(self):
        try:
            return float(self.config_parser.get('SIMULATION', 'FRICTION_TRANSITION_REYNOLDS'))
        except:
            return 1.0e6

    def getFrictionMaxCf(self):
        try:
            return float(self.config_parser.get('SIMULATION', 'FRICTION_MAX_CF'))
        except:
            return 0.1

    def getFrictionDriverCf(self):
        try:
            return float(self.config_parser.get('SIMULATION', 'FRICTION_DRIVER_CF'))
        except:
            return self.getFrictionCoefficient()

    def getShockDetectionThreshold(self):
        try:
            return float(self.config_parser.get('SIMULATION', 'SHOCK_DETECTION_THRESHOLD'))
        except:
            return 0.05 
    
    def isWallHeatTransferActive(self):
        try:
            res = str(self.config_parser.get('SIMULATION', 'WALL_HEAT_TRANSFER_ACTIVE')).lower()
            if res=='yes' or res=='true':
                return True
            else:
                return False
        except:
            return False 
    
    def getWallHeatFlux(self):
        try:
            return float(self.config_parser.get('SIMULATION', 'WALL_HEAT_FLUX'))
        except:
            return 0.0 # default value

    def getNozzleFilePath(self):
        return str(self.config_parser.get('GEOMETRY', 'NOZZLE_FILEPATH'))
    
    
    def getAreaReference(self):
        try:
            return float(self.config_parser.get('GEOMETRY', 'REFERENCE_AREA')) 
        except:
            return 1.0 # default
    
    
    def getOutputFrequency(self):
        try:
            return int(self.config_parser.get('OUTPUT', 'OUTPUT_FREQUENCY')) 
        except:
            return 250 # default value

    def isLookupTableActive(self):
        """Check whether 2D Look-Up Table (LuT) acceleration is enabled."""
        for section in ('SIMULATION', 'FLUID'):
            try:
                res = str(self.config_parser.get(section, 'USE_LUT')).lower()
                return res in ['yes', 'true', '1']
            except Exception:
                pass
        return False

    def getLookupTableGridSize(self):
        """Return 2D grid resolution (nP, nT) for Look-Up Table (default: (250, 250), max: (1000, 1000))."""
        for section in ('FLUID', 'SIMULATION'):
            try:
                res = str(self.config_parser.get(section, 'LUT_GRID_SIZE')).strip()
                parts = [int(v.strip()) for v in res.split(',')]
                if len(parts) == 1:
                    val = min(1000, max(20, parts[0]))
                    return (val, val)
                nP = min(1000, max(20, parts[0]))
                nT = min(1000, max(20, parts[1]))
                return (nP, nT)
            except Exception:
                pass
        return (250, 250)

    def getLookupTablePressureRange(self):
        """Return (P_min, P_max) in Pascals for Look-Up Table domain."""
        for section in ('FLUID', 'SIMULATION'):
            try:
                p_min = float(self.config_parser.get(section, 'LUT_PRESSURE_MIN'))
                p_max = float(self.config_parser.get(section, 'LUT_PRESSURE_MAX'))
                return p_min, p_max
            except Exception:
                pass
        # Default estimation from initial states
        pL = self.getPressureLeft()
        pR = self.getPressureRight()
        return min(pL, pR) * 0.5, max(pL, pR) * 1.5

    def getLookupTableTemperatureRange(self):
        """Return (T_min, T_max) in Kelvin for Look-Up Table domain."""
        for section in ('FLUID', 'SIMULATION'):
            try:
                T_min = float(self.config_parser.get(section, 'LUT_TEMPERATURE_MIN'))
                T_max = float(self.config_parser.get(section, 'LUT_TEMPERATURE_MAX'))
                return T_min, T_max
            except Exception:
                pass
        # Default estimation for shock tube real fluids
        try:
            TL = self.getTemperatureLeft()
            TR = self.getTemperatureRight()
            return min(TL, TR) * 0.7, max(TL, TR) * 1.3
        except Exception:
            return 200.0, 2000.0