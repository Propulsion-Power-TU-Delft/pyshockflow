import CoolProp.CoolProp as CP
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import sys
import fluid_properties.fluid_properties as FP

class FluidIdeal():
    """
    Ideal Fluid Class, where thermodynamic properties and transformation are computed with ideal gas laws
    """
    def __init__(self, gmma, Rgas):
        self.gmma = gmma
        self.Rgas = Rgas
    
    def computeStaticEnergy_p_rho(self, p, rho):
        return (p / (self.gmma - 1) / rho)
    
    def computePressure_rho_e(self, rho, e):
        return (self.gmma-1)*rho*e
    
    def computeSoundSpeed_p_rho(self, p, rho):
        return np.sqrt(self.gmma*p/rho)
    
    def computeMach_u_p_rho(self, u, p, rho):
        soundSpeed = self.computeSoundSpeed_p_rho(p, rho)
        return np.abs(u)/soundSpeed

    def computeTemperature_p_rho(self, p, rho):
        return (p/rho)/self.Rgas

    def computeDensity_p_T(self, p, T):
        return p/self.Rgas/T

    def computeEntropy_p_rho(self, p, rho):
        return p/(rho**self.gmma)

    def computeFunDerGamma_p_rho(self, p, rho):
        if isinstance(p, np.ndarray): # handle the case when the inputs are arrays
            return 0.5*(self.gmma+1)+np.zeros_like(p)
        else:
            return 0.5*(self.gmma+1)

    def computeComprFactorZ_p_rho(self, p, rho):
        if isinstance(p, np.ndarray):
            return 1+np.zeros_like(p)
        else:
            return 1

    def computeTotalPressure_p_M(self, p, M):
        return p*(1+(self.gmma-1)/2*M**2)**(self.gmma/(self.gmma-1))

    def computeMach_pt_p(self, pt, p):
        mach = np.sqrt( 2/(self.gmma-1) * ((pt/p)**((self.gmma-1)/self.gmma)-1) )
        return mach

    def computeTotalTemperature_T_M(self, T, M):
        return T*(1+(self.gmma-1)/2*M**2)

    def computeTemperature_Tt_M(self, Tt, M):
        return Tt/(1+(self.gmma-1)/2*M**2)

    def computePressure_Pt_M(self, Pt, M):
        return Pt/((1+(self.gmma-1)/2*M**2)**(self.gmma/(self.gmma-1)))
    
    def computeInletQuantities(self, pressure, totPressure, totTemperature, direction):
        mach = self.computeMach_pt_p(totPressure, pressure)
        temperature = self.computeTemperature_Tt_M(totTemperature, mach)
        density = self.computeDensity_p_T(pressure, temperature)
        soundSpeed = self.computeSoundSpeed_p_rho(pressure, density)
        velocity = mach*soundSpeed*direction
        energy = self.computeStaticEnergy_p_rho(pressure, density)
        return density, velocity, energy

    def compute_gammapv_p_rho(self, p, rho):
        if isinstance(p, np.ndarray):
            gmma_pv = np.zeros_like(p)+self.gmma
        else:
            gmma_pv = self.gmma
        return gmma_pv

    def computeChiKappa_VinokurScheme_p_rho(self, p, rho):
        chi = 0
        kappa = self.gmma-1
        return chi, kappa

    def compute_e_and_a_p_rho(self, p, rho):
        e = p / ((self.gmma - 1.0) * rho)
        a = np.sqrt(self.gmma * p / rho)
        return e, a

    def compute_e_and_chi_kappa_p_rho(self, p, rho):
        e = p / ((self.gmma - 1.0) * rho)
        gm1 = self.gmma - 1.0
        if isinstance(p, np.ndarray):
            chi = np.zeros_like(p)
            kappa = np.full_like(p, gm1)
        else:
            chi = 0.0
            kappa = gm1
        return e, chi, kappa
    
    
class FluidReal():
    """
    Real Fluid Class, where thermodynamic properties and transformations are taken from CoolProp / external backends
    """
    def __init__(self, fluid_name, fluid_library, print_error=True):
        self.fluid_name = fluid_name
        self.fluid_library = fluid_library
        self.fluid = FP.fluid(fluid_library, fluid_name, print_error=print_error)
        
        # Initialize persistent low-level AbstractState if available for direct high-speed evaluation
        self._backend = None
        if fluid_library in ['CoolProp', 'HEOS']:
            try:
                import CoolProp
                from CoolProp.CoolProp import AbstractState
                self._backend = AbstractState('HEOS', fluid_name)
            except Exception:
                self._backend = None

    def __getstate__(self):
        state = self.__dict__.copy()
        state['_backend'] = None
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        if self.fluid_library in ['CoolProp', 'HEOS']:
            try:
                import CoolProp
                from CoolProp.CoolProp import AbstractState
                self._backend = AbstractState('HEOS', self.fluid_name)
            except Exception:
                self._backend = None
        else:
            self._backend = None

    def computeStaticEnergy_p_rho(self, p, rho):
        if np.isscalar(p):
            if self._backend is not None:
                try:
                    import CoolProp
                    self._backend.update(CoolProp.DmassP_INPUTS, float(rho), float(p))
                    return self._backend.umass()
                except Exception:
                    pass
            return FP.PropsSI('U', 'P', p, 'D', rho, self.fluid)
        # Array inputs
        if self._backend is not None:
            import CoolProp
            out = np.empty(len(p))
            for i in range(len(p)):
                try:
                    self._backend.update(CoolProp.DmassP_INPUTS, float(rho[i]), float(p[i]))
                    out[i] = self._backend.umass()
                except Exception:
                    out[i] = FP.PropsSI('U', 'P', p[i], 'D', rho[i], self.fluid)
            return out
        return FP.PropsSI('U', 'P', p, 'D', rho, self.fluid)
    
    def computePressure_rho_e(self, rho, e):
        if np.isscalar(rho):
            if self._backend is not None:
                try:
                    import CoolProp
                    self._backend.update(CoolProp.DmassUmass_INPUTS, float(rho), float(e))
                    return self._backend.p()
                except Exception:
                    pass
            return FP.PropsSI('P', 'D', rho, 'U', e, self.fluid)
        # Array inputs
        if self._backend is not None:
            import CoolProp
            out = np.empty(len(rho))
            for i in range(len(rho)):
                try:
                    self._backend.update(CoolProp.DmassUmass_INPUTS, float(rho[i]), float(e[i]))
                    out[i] = self._backend.p()
                except Exception:
                    out[i] = FP.PropsSI('P', 'D', rho[i], 'U', e[i], self.fluid)
            return out
        return FP.PropsSI('P', 'D', rho, 'U', e, self.fluid)

    def computeSoundSpeed_p_rho(self, p, rho):
        if np.isscalar(p):
            if self._backend is not None:
                try:
                    import CoolProp
                    self._backend.update(CoolProp.DmassP_INPUTS, float(rho), float(p))
                    return self._backend.speed_sound()
                except Exception:
                    pass
            try:
                a = FP.PropsSI("A", "P", p, "D", rho, self.fluid)
                return a
            except Exception:
                return self._computeSoundSpeed_twophase(p, rho)
        
        # Array inputs
        out = np.empty(len(p))
        for i in range(len(p)):
            out[i] = self.computeSoundSpeed_p_rho(p[i], rho[i])
        return out

    def _computeSoundSpeed_twophase(self, p, rho):
        T = self.computeTemperature_p_rho(p, rho)
        try:
            Q = FP.PropsSI("Q", "T", T, "P", p, self.fluid)
        except Exception:
            Q = 1
        a_liquid = FP.PropsSI("A", "T", T, "Q", 0, self.fluid)
        a_vapor = FP.PropsSI("A", "T", T, "Q", 1, self.fluid)
        return (1 - Q) * a_liquid + Q * a_vapor

    def compute_e_and_a_p_rho(self, p, rho):
        """Simultaneously evaluate internal energy e and sound speed a from (p, rho).
        Performs a single state update rather than two separate thermodynamic solves.
        """
        if np.isscalar(p):
            if self._backend is not None:
                try:
                    import CoolProp
                    self._backend.update(CoolProp.DmassP_INPUTS, float(rho), float(p))
                    return self._backend.umass(), self._backend.speed_sound()
                except Exception:
                    pass
            return self.computeStaticEnergy_p_rho(p, rho), self.computeSoundSpeed_p_rho(p, rho)

        n = len(p)
        e = np.empty(n)
        a = np.empty(n)
        if self._backend is not None:
            import CoolProp
            for i in range(n):
                try:
                    self._backend.update(CoolProp.DmassP_INPUTS, float(rho[i]), float(p[i]))
                    e[i] = self._backend.umass()
                    a[i] = self._backend.speed_sound()
                except Exception:
                    e[i] = self.computeStaticEnergy_p_rho(p[i], rho[i])
                    a[i] = self.computeSoundSpeed_p_rho(p[i], rho[i])
        else:
            for i in range(n):
                e[i] = self.computeStaticEnergy_p_rho(p[i], rho[i])
                a[i] = self.computeSoundSpeed_p_rho(p[i], rho[i])
        return e, a

    def compute_e_and_chi_kappa_p_rho(self, p, rho):
        """Simultaneously evaluate internal energy e, and Vinokur derivatives (chi, kappa).
        Uses exact analytical partial derivatives from AbstractState when available.
        """
        if np.isscalar(p):
            if self._backend is not None:
                try:
                    import CoolProp
                    self._backend.update(CoolProp.DmassP_INPUTS, float(rho), float(p))
                    e = self._backend.umass()
                    dp_drho_econst = self._backend.first_partial_deriv(CoolProp.iP, CoolProp.iDmass, CoolProp.iUmass)
                    dp_de_rhoconst = self._backend.first_partial_deriv(CoolProp.iP, CoolProp.iUmass, CoolProp.iDmass)
                    chi = dp_drho_econst - (e / rho) * dp_de_rhoconst
                    kappa = dp_de_rhoconst / rho
                    return e, chi, kappa
                except Exception:
                    pass
            e = FP.PropsSI("U", "P", p, "D", rho, self.fluid)
            dp_drho_econst = FP.PropsSI("d(P)/d(D)|U", "P", p, "D", rho, self.fluid)
            dp_de_rhoconst = FP.PropsSI("d(P)/d(U)|D", "P", p, "D", rho, self.fluid)
            chi = dp_drho_econst - e/rho * dp_de_rhoconst
            kappa = dp_de_rhoconst / rho
            return e, chi, kappa

        n = len(p)
        e = np.empty(n)
        chi = np.empty(n)
        kappa = np.empty(n)
        if self._backend is not None:
            import CoolProp
            for i in range(n):
                try:
                    self._backend.update(CoolProp.DmassP_INPUTS, float(rho[i]), float(p[i]))
                    e[i] = self._backend.umass()
                    dp_drho_econst = self._backend.first_partial_deriv(CoolProp.iP, CoolProp.iDmass, CoolProp.iUmass)
                    dp_de_rhoconst = self._backend.first_partial_deriv(CoolProp.iP, CoolProp.iUmass, CoolProp.iDmass)
                    chi[i] = dp_drho_econst - (e[i] / rho[i]) * dp_de_rhoconst
                    kappa[i] = dp_de_rhoconst / rho[i]
                except Exception:
                    e[i] = FP.PropsSI("U", "P", p[i], "D", rho[i], self.fluid)
                    dp_drho_econst = FP.PropsSI("d(P)/d(D)|U", "P", p[i], "D", rho[i], self.fluid)
                    dp_de_rhoconst = FP.PropsSI("d(P)/d(U)|D", "P", p[i], "D", rho[i], self.fluid)
                    chi[i] = dp_drho_econst - e[i]/rho[i] * dp_de_rhoconst
                    kappa[i] = dp_de_rhoconst / rho[i]
        else:
            for i in range(n):
                e[i] = FP.PropsSI("U", "P", p[i], "D", rho[i], self.fluid)
                dp_drho_econst = FP.PropsSI("d(P)/d(D)|U", "P", p[i], "D", rho[i], self.fluid)
                dp_de_rhoconst = FP.PropsSI("d(P)/d(U)|D", "P", p[i], "D", rho[i], self.fluid)
                chi[i] = dp_drho_econst - e[i]/rho[i] * dp_de_rhoconst
                kappa[i] = dp_de_rhoconst / rho[i]
        return e, chi, kappa

    def computeMach_u_p_rho(self, u, p, rho):
        soundSpeed = self.computeSoundSpeed_p_rho(p, rho)
        return np.abs(u)/soundSpeed

    def computeTemperature_p_rho(self, p, rho):
        T = FP.PropsSI('T', 'P', p, 'D', rho, self.fluid)
        return T

    def computeDensity_p_T(self, p, T):
        rho = FP.PropsSI('D', 'P', p, 'T', T, self.fluid)
        return rho

    def computeEntropy_p_rho(self, p, rho):
        s = FP.PropsSI('S', 'P', p, 'D', rho, self.fluid)
        return s

    def computeEntropy_p_T(self, p, T):
        s = FP.PropsSI('S', 'P', p, 'T', T, self.fluid)
        return s

    def computeFunDerGamma_p_rho(self, p, rho):
        try: # if single phase this will work
            G = FP.PropsSI("FUNDAMENTAL_DERIVATIVE_OF_GAS_DYNAMICS", "P", p, "D", rho, self.fluid)
            return G
        except: # if close to two phase, we need to do like the speed of sound
            T = self.computeTemperature_p_rho(p, rho)
            try:
                Q = FP.PropsSI("Q", "T", T, "P", p, self.fluid)
            except:
                # if the state is very close to saturation line it fails to find the quality -> set artifically to 1
                Q = 1

            # G in liquid and vapor phases at the given T
            G_liquid = FP.PropsSI("FUNDAMENTAL_DERIVATIVE_OF_GAS_DYNAMICS", "T", T, "Q", 0, self.fluid)  # sound speed for liquid phase
            G_vapor = FP.PropsSI("FUNDAMENTAL_DERIVATIVE_OF_GAS_DYNAMICS", "T", T, "Q", 1, self.fluid)   # sound speed for vapor phase

            # Calculate weighted G based on quality
            G = (1 - Q) * G_liquid + Q * G_vapor
            return G

    def computeComprFactorZ_p_rho(self, p, rho):
        Z = FP.PropsSI('Z', 'P', p, 'D', rho, self.fluid)
        return Z

    
    def computeInletQuantities(self, pressure, totPressure, totTemperature, direction):
        """The full state must be reconstructed from the quantities given in the arguments.
        The entropy of the static and total state must be the same by definition. This is used to find the temperature.

        Args:
            pressure (float): static pressure
            totPressure (float): total pressure
            totTemperature (float): total temperature
        """
        def compute_function_residual(temperatureGuess):
            entropyStatic = self.computeEntropy_p_T(pressure, temperatureGuess)
            entropyTotal = self.computeEntropy_p_T(totPressure, totTemperature)
            residual = entropyStatic - entropyTotal
            return residual

        # temperature = fsolve(compute_function_residual, totTemperature, xtol=1e-8)[0]
        temperature, info, ier, msg = fsolve(
            compute_function_residual,
            totTemperature,
            xtol=1e-6,
            full_output=True
        )
        if ier != 1:
            raise RuntimeError(f"fsolve did not converge: {msg}")
        
        temperature = temperature[0]
        density = self.computeDensity_p_T(pressure, temperature)
        gamma_pv = self.compute_gammapv_p_rho(pressure, density)
        mach = self.computeMach_pt_p_gammapv(totPressure, pressure, gamma_pv)
        soundSpeed = self.computeSoundSpeed_p_rho(pressure, density)
        velocity = direction * mach * soundSpeed
        energy = self.computeStaticEnergy_p_rho(pressure, density)
        return density, velocity, energy


    def compute_gammapv_p_rho(self, p, rho):
        cp = FP.PropsSI("Cpmass", "P", p, "D", rho, self.fluid)
        cv = FP.PropsSI("Cvmass", "P", p, "D", rho, self.fluid)
        dp_drho_T = FP.PropsSI("d(P)/d(D)|T", "P", p, "D", rho, self.fluid)
        dp_dv_T = - rho**2 * dp_drho_T
        gmma_pv = -1/(p*rho) * cp/cv * dp_dv_T
        return gmma_pv


    def compute_gammapt_p_T(self, p, T):
        rho = FP.PropsSI("D", "P", p, "T", T, self.fluid)
        d_rho_dT_P = FP.PropsSI("d(D)/d(T)|P", "P", p, "T", T, self.fluid)
        dv_dT_P = - d_rho_dT_P / (rho**2)
        cp = FP.PropsSI("Cpmass", "P", p, "T", T, self.fluid)
        gamma_pT = 1 / (1 - p/cp*dv_dT_P)
        return gamma_pT


    def computeMach_pt_p_gammapv(self, pt, p, gamma_pv):
        """Reference to equation 8.10 Nederstigt MS thesis
        """
        mach = np.sqrt(2/(gamma_pv-1) * ((pt/p)**((gamma_pv-1)/gamma_pv) - 1))
        return mach
    

    def computeChiKappa_VinokurScheme_p_rho(self, p, rho):
        e = FP.PropsSI("U", "P", p, "D", rho, self.fluid)
        dp_drho_econst = FP.PropsSI("d(P)/d(D)|U", "P", p, "D", rho, self.fluid)
        dp_de_rhoconst = FP.PropsSI("d(P)/d(U)|D", "P", p, "D", rho, self.fluid)
        chi = dp_drho_econst - e/rho * dp_de_rhoconst
        kappa = dp_de_rhoconst / rho
        return chi, kappa
        
        
        

            