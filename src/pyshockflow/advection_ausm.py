import numpy as np
from pyshockflow import FluidIdeal


class AdvectionAUSMplusUP:
    """
    AUSM+-up (Advection Upstream Splitting Method for all speeds) Riemann solver
    for the 1D Euler equations.
    
    Supports both ideal gases and general real gas equations of state (FluidIdeal and FluidReal).
    
    References:
    - Liou, M.-S., "A sequel to AUSM, Part II: AUSM+-up for all speeds",
      Journal of Computational Physics, 214(1), 137-170, 2006.
    - Edwards, J. R., & Liou, M.-S., "Low-Diffusion Flux-Splitting Methods for Real Fluid Flows",
      AIAA Journal, 1998.
    """
    def __init__(self, rhoL, rhoR, uL, uR, pL, pR, fluid,
                 Kp=0.25, Ku=0.75, sigma=1.0, beta=0.125):
        self.rhoL = float(rhoL)
        self.rhoR = float(rhoR)
        self.uL = float(uL)
        self.uR = float(uR)
        self.pL = float(pL)
        self.pR = float(pR)
        self.fluid = fluid

        self.Kp = Kp
        self.Ku = Ku
        self.sigma = sigma
        self.beta = beta

        self.eL = fluid.computeStaticEnergy_p_rho(self.pL, self.rhoL)
        self.eR = fluid.computeStaticEnergy_p_rho(self.pR, self.rhoR)
        self.aL = fluid.computeSoundSpeed_p_rho(self.pL, self.rhoL)
        self.aR = fluid.computeSoundSpeed_p_rho(self.pR, self.rhoR)

        # Total specific enthalpies: H = e + 0.5*u^2 + p/rho
        self.HL = self.eL + 0.5 * self.uL**2 + self.pL / self.rhoL
        self.HR = self.eR + 0.5 * self.uR**2 + self.pR / self.rhoR

    def _split_mach_plus(self, M):
        if abs(M) <= 1.0:
            return 0.25 * (M + 1.0)**2 + self.beta * (M**2 - 1.0)**2
        return 0.5 * (M + abs(M))

    def _split_mach_minus(self, M):
        if abs(M) <= 1.0:
            return -0.25 * (M - 1.0)**2 - self.beta * (M**2 - 1.0)**2
        return 0.5 * (M - abs(M))

    def _split_pressure_plus(self, M, alpha):
        if abs(M) <= 1.0:
            return 0.25 * (M + 1.0)**2 * (2.0 - M) + alpha * M * (M**2 - 1.0)**2
        return 0.5 * (1.0 + np.sign(M))

    def _split_pressure_minus(self, M, alpha):
        if abs(M) <= 1.0:
            return 0.25 * (M - 1.0)**2 * (2.0 + M) - alpha * M * (M**2 - 1.0)**2
        return 0.5 * (1.0 - np.sign(M))

    def computeFlux(self, entropyFixActive=False, fixCoefficient=0.2):
        """
        Compute the AUSM+-up interface flux vector.
        Arguments entropyFixActive and fixCoefficient are accepted for interface compatibility.
        """
        # Interface sound speed (arithmetic mean for general EOS)
        a_half = 0.5 * (self.aL + self.aR)

        # Left / right Mach numbers
        ML = self.uL / a_half
        MR = self.uR / a_half

        # All-speed reference Mach scaling (accounts for both convective velocity and pressure shocks)
        M_bar_sq = (self.uL**2 + self.uR**2) / (2.0 * a_half**2)
        dp_ratio = abs(self.pR - self.pL) / (self.pL + self.pR + 1e-30)
        M0_sq = min(1.0, max(M_bar_sq, dp_ratio**2, 1e-2))
        M0 = np.sqrt(M0_sq)
        fa = max(M0 * (2.0 - M0), 0.1)

        alpha = (3.0 / 16.0) * (-4.0 + 5.0 * fa**2)

        # Mach and Pressure polynomial splittings
        M_plus = self._split_mach_plus(ML)
        M_minus = self._split_mach_minus(MR)
        P_plus = self._split_pressure_plus(ML, alpha)
        P_minus = self._split_pressure_minus(MR, alpha)

        # Pressure diffusion term for interface Mach number (-up term)
        rho_half = 0.5 * (self.rhoL + self.rhoR)
        Dp = (self.Kp / fa) * max(1.0 - self.sigma * M_bar_sq, 0.0) * (self.pR - self.pL) / (rho_half * a_half**2)
        M_half = M_plus + M_minus - Dp

        # Velocity diffusion term for interface pressure (-up term)
        Du = self.Ku * P_plus * P_minus * (self.rhoL + self.rhoR) * fa * a_half * (self.uR - self.uL)
        p_half = P_plus * self.pL + P_minus * self.pR - Du

        # Interface mass flux
        if M_half >= 0.0:
            m_dot = a_half * M_half * self.rhoL
            f1 = m_dot
            f2 = m_dot * self.uL + p_half
            f3 = m_dot * self.HL
        else:
            m_dot = a_half * M_half * self.rhoR
            f1 = m_dot
            f2 = m_dot * self.uR + p_half
            f3 = m_dot * self.HR

        return np.array([f1, f2, f3])
