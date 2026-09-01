"""
Vectorized AUSM+-up approximate Riemann solvers for the 1D Euler equations.

These functions compute numerical fluxes at ALL interfaces simultaneously
using NumPy array operations.

Two variants are provided:
- ``compute_ausm_flux_ideal``: AUSM+-up for ideal gas
- ``compute_ausm_flux_real``: AUSM+-up for real gases (general EOS)
"""

import numpy as np


def compute_ausm_flux_ideal(rhoL, rhoR, uL, uR, pL, pR, gmma,
                            Kp=0.25, Ku=0.75, sigma=1.0, beta=0.125):
    """Vectorized AUSM+-up flux for ideal gas.

    Parameters
    ----------
    rhoL, rhoR, uL, uR, pL, pR : ndarray (N,)
        Left / right primitive states at each of the *N* interfaces.
    gmma : float
        Ratio of specific heats γ.
    Kp, Ku, sigma, beta : float, optional

    Returns
    -------
    flux : ndarray (N, 3)
        Numerical flux at each interface.
    """
    gm1 = gmma - 1.0

    # Static energy (ideal gas)
    eL = pL / (gm1 * rhoL)
    eR = pR / (gm1 * rhoR)

    # Sound speed (ideal gas)
    aL = np.sqrt(gmma * pL / rhoL)
    aR = np.sqrt(gmma * pR / rhoR)

    return compute_ausm_flux_real(rhoL, rhoR, uL, uR, pL, pR, eL, eR, aL, aR,
                                  Kp=Kp, Ku=Ku, sigma=sigma, beta=beta)


def compute_ausm_flux_real(rhoL, rhoR, uL, uR, pL, pR, eL, eR, aL, aR,
                           Kp=0.25, Ku=0.75, sigma=1.0, beta=0.125):
    """Vectorized AUSM+-up flux for real gas (general equation of state).

    Parameters
    ----------
    rhoL, rhoR, uL, uR, pL, pR : ndarray (N,)
        Left / right primitive states at each of the *N* interfaces.
    eL, eR : ndarray (N,)
        Static internal energy per unit mass at left / right states.
    aL, aR : ndarray (N,)
        Sound speed at left / right states.
    Kp, Ku, sigma, beta : float, optional

    Returns
    -------
    flux : ndarray (N, 3)
        Numerical flux at each interface.
    """
    # Total specific enthalpies: H = e + 0.5*u^2 + p/rho
    HL = eL + 0.5 * uL**2 + pL / rhoL
    HR = eR + 0.5 * uR**2 + pR / rhoR

    # Interface sound speed (arithmetic mean for general EOS)
    a_half = 0.5 * (aL + aR)

    # Mach numbers
    ML = uL / a_half
    MR = uR / a_half

    # All-speed reference Mach scaling
    M_bar_sq = (uL**2 + uR**2) / (2.0 * a_half**2)
    dp_ratio = np.abs(pR - pL) / (pL + pR + 1e-30)
    M0_sq = np.clip(np.maximum(M_bar_sq, dp_ratio**2), 1e-2, 1.0)
    M0 = np.sqrt(M0_sq)
    fa = np.maximum(M0 * (2.0 - M0), 0.1)

    alpha = (3.0 / 16.0) * (-4.0 + 5.0 * fa**2)

    # Mach polynomial splitting
    abs_ML = np.abs(ML)
    M_plus = np.where(abs_ML <= 1.0,
                      0.25 * (ML + 1.0)**2 + beta * (ML**2 - 1.0)**2,
                      0.5 * (ML + abs_ML))

    abs_MR = np.abs(MR)
    M_minus = np.where(abs_MR <= 1.0,
                       -0.25 * (MR - 1.0)**2 - beta * (MR**2 - 1.0)**2,
                       0.5 * (MR - abs_MR))

    # Pressure polynomial splitting
    P_plus = np.where(abs_ML <= 1.0,
                      0.25 * (ML + 1.0)**2 * (2.0 - ML) + alpha * ML * (ML**2 - 1.0)**2,
                      0.5 * (1.0 + np.sign(ML)))

    P_minus = np.where(abs_MR <= 1.0,
                       0.25 * (MR - 1.0)**2 * (2.0 + MR) - alpha * MR * (MR**2 - 1.0)**2,
                       0.5 * (1.0 - np.sign(MR)))

    # Pressure diffusion term for interface Mach number (-up term)
    rho_half = 0.5 * (rhoL + rhoR)
    Dp = (Kp / fa) * np.maximum(1.0 - sigma * M_bar_sq, 0.0) * (pR - pL) / (rho_half * a_half**2)
    M_half = M_plus + M_minus - Dp

    # Velocity diffusion term for interface pressure (-up term)
    Du = Ku * P_plus * P_minus * (rhoL + rhoR) * fa * a_half * (uR - uL)
    p_half = P_plus * pL + P_minus * pR - Du

    # Numerical flux assembly
    m_dot = np.where(M_half >= 0.0, a_half * M_half * rhoL, a_half * M_half * rhoR)
    f1 = m_dot
    f2 = np.where(M_half >= 0.0, m_dot * uL + p_half, m_dot * uR + p_half)
    f3 = np.where(M_half >= 0.0, m_dot * HL, m_dot * HR)

    return np.column_stack((f1, f2, f3))
