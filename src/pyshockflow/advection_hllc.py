import numpy as np
from pyshockflow import FluidIdeal
from pyshockflow.math_utils import getConservativesFromPrimitives


class AdvectionHLLC:
    """
    HLLC (Harten-Lax-van Leer Contact) approximate Riemann solver for 1D Euler equations.
    Supports both ideal gases and general real gas equations of state (FluidIdeal and FluidReal).
    
    References:
    - Toro, E. F., "Riemann Solvers and Numerical Methods for Fluid Dynamics", Springer, 2009.
    - Batten, P., Clarke, N., Lambert, C., & Causon, D. M., "On the Choice of Wave Speeds for the HLLC
      Riemann Solver", SIAM J. Sci. Comput., 1997.
    """
    def __init__(self, rhoL, rhoR, uL, uR, pL, pR, fluid):
        self.rhoL = float(rhoL)
        self.rhoR = float(rhoR)
        self.uL = float(uL)
        self.uR = float(uR)
        self.pL = float(pL)
        self.pR = float(pR)
        self.fluid = fluid

        self.eL = fluid.computeStaticEnergy_p_rho(self.pL, self.rhoL)
        self.eR = fluid.computeStaticEnergy_p_rho(self.pR, self.rhoR)
        self.aL = fluid.computeSoundSpeed_p_rho(self.pL, self.rhoL)
        self.aR = fluid.computeSoundSpeed_p_rho(self.pR, self.rhoR)

        # Conservative state vectors: U = [rho, rho*u, rho*E]
        EL = self.eL + 0.5 * self.uL**2
        ER = self.eR + 0.5 * self.uR**2
        self.UL = np.array([self.rhoL, self.rhoL * self.uL, self.rhoL * EL])
        self.UR = np.array([self.rhoR, self.rhoR * self.uR, self.rhoR * ER])

        # Physical Euler fluxes: F = [rho*u, rho*u^2 + p, u*(rho*E + p)]
        self.FL = np.array([
            self.rhoL * self.uL,
            self.rhoL * self.uL**2 + self.pL,
            self.uL * (self.rhoL * EL + self.pL)
        ])
        self.FR = np.array([
            self.rhoR * self.uR,
            self.rhoR * self.uR**2 + self.pR,
            self.uR * (self.rhoR * ER + self.pR)
        ])

    def computeWaveSpeeds(self):
        """
        Estimate left wave speed SL, right wave speed SR, and contact wave speed S_star.
        Uses Davis / Einfeldt wave speed estimates suitable for general EOS.
        """
        # Davis wave speed estimates
        self.SL = min(self.uL - self.aL, self.uR - self.aR)
        self.SR = max(self.uL + self.aL, self.uR + self.aR)

        # Intermediate contact wave speed S*
        num = (self.pR - self.pL
               + self.rhoL * self.uL * (self.SL - self.uL)
               - self.rhoR * self.uR * (self.SR - self.uR))
        den = self.rhoL * (self.SL - self.uL) - self.rhoR * (self.SR - self.uR)

        if abs(den) < 1e-14:
            self.S_star = 0.5 * (self.uL + self.uR)
        else:
            self.S_star = num / den

    def computeStarState(self, side='L'):
        """
        Compute star state U_*L or U_*R for the HLLC Riemann solver.
        """
        if side == 'L':
            rho, u, p, e, S, U = self.rhoL, self.uL, self.pL, self.eL, self.SL, self.UL
        else:
            rho, u, p, e, S, U = self.rhoR, self.uR, self.pR, self.eR, self.SR, self.UR

        E = e + 0.5 * u**2
        factor = rho * (S - u) / (S - self.S_star + 1e-30)

        # Star energy per unit mass: E_* = E + (S* - u) * [S* + p / (rho * (S - u))]
        E_star = E + (self.S_star - u) * (self.S_star + p / (rho * (S - u) + 1e-30))

        U_star = np.array([
            factor,
            factor * self.S_star,
            factor * E_star
        ])
        return U_star

    def computeFlux(self, entropyFixActive=False, fixCoefficient=0.2):
        """
        Compute the numerical interface flux using the HLLC formulation.
        Note: HLLC naturally satisfies the entropy condition and does not require an entropy fix.
        Arguments entropyFixActive and fixCoefficient are accepted for interface compatibility.
        """
        self.computeWaveSpeeds()

        if self.SL >= 0.0:
            return self.FL.copy()
        elif self.S_star >= 0.0:
            U_starL = self.computeStarState('L')
            return self.FL + self.SL * (U_starL - self.UL)
        elif self.SR >= 0.0:
            U_starR = self.computeStarState('R')
            return self.FR + self.SR * (U_starR - self.UR)
        else:
            return self.FR.copy()
