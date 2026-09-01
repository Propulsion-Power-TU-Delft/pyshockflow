"""
Vectorized HLLC approximate Riemann solvers for the 1D Euler equations.

These functions compute numerical fluxes at ALL interfaces simultaneously
using NumPy array operations.

Two variants are provided:
- ``compute_hllc_flux_ideal``: HLLC for ideal gas
- ``compute_hllc_flux_real``: HLLC for real gases (general EOS)
"""

import numpy as np


def compute_hllc_flux_ideal(rhoL, rhoR, uL, uR, pL, pR, gmma):
    """Vectorized HLLC flux for ideal gas.

    Parameters
    ----------
    rhoL, rhoR, uL, uR, pL, pR : ndarray (N,)
        Left / right primitive states at each of the *N* interfaces.
    gmma : float
        Ratio of specific heats γ.

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

    return compute_hllc_flux_real(rhoL, rhoR, uL, uR, pL, pR, eL, eR, aL, aR)


def compute_hllc_flux_real(rhoL, rhoR, uL, uR, pL, pR, eL, eR, aL, aR):
    """Vectorized HLLC flux for real gas (general equation of state).

    Parameters
    ----------
    rhoL, rhoR, uL, uR, pL, pR : ndarray (N,)
        Left / right primitive states at each of the *N* interfaces.
    eL, eR : ndarray (N,)
        Static internal energy per unit mass at left / right states.
    aL, aR : ndarray (N,)
        Sound speed at left / right states.

    Returns
    -------
    flux : ndarray (N, 3)
        Numerical flux at each interface.
    """
    # Total energy per unit mass
    EL = eL + 0.5 * uL**2
    ER = eR + 0.5 * uR**2

    # Physical Euler fluxes
    FL1 = rhoL * uL
    FL2 = rhoL * uL**2 + pL
    FL3 = uL * (rhoL * EL + pL)

    FR1 = rhoR * uR
    FR2 = rhoR * uR**2 + pR
    FR3 = uR * (rhoR * ER + pR)

    # Davis wave speed estimates
    SL = np.minimum(uL - aL, uR - aR)
    SR = np.maximum(uL + aL, uR + aR)

    # Contact wave speed S*
    num = pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR)
    den = rhoL * (SL - uL) - rhoR * (SR - uR)

    # Safe division for contact wave speed
    safe_den = np.where(np.abs(den) < 1e-14, 1.0, den)
    S_star = np.where(np.abs(den) < 1e-14, 0.5 * (uL + uR), num / safe_den)

    # Left star state
    diff_L = SL - S_star
    safe_diff_L = np.where(np.abs(diff_L) < 1e-30, 1e-30, diff_L)
    diff_SL_uL = SL - uL
    safe_SL_uL = np.where(np.abs(diff_SL_uL) < 1e-30, 1e-30, diff_SL_uL)

    factor_L = rhoL * diff_SL_uL / safe_diff_L
    E_star_L = EL + (S_star - uL) * (S_star + pL / (rhoL * safe_SL_uL))

    U_star_L1 = factor_L
    U_star_L2 = factor_L * S_star
    U_star_L3 = factor_L * E_star_L

    F_star_L1 = FL1 + SL * (U_star_L1 - rhoL)
    F_star_L2 = FL2 + SL * (U_star_L2 - rhoL * uL)
    F_star_L3 = FL3 + SL * (U_star_L3 - rhoL * EL)

    # Right star state
    diff_R = SR - S_star
    safe_diff_R = np.where(np.abs(diff_R) < 1e-30, 1e-30, diff_R)
    diff_SR_uR = SR - uR
    safe_SR_uR = np.where(np.abs(diff_SR_uR) < 1e-30, 1e-30, diff_SR_uR)

    factor_R = rhoR * diff_SR_uR / safe_diff_R
    E_star_R = ER + (S_star - uR) * (S_star + pR / (rhoR * safe_SR_uR))

    U_star_R1 = factor_R
    U_star_R2 = factor_R * S_star
    U_star_R3 = factor_R * E_star_R

    F_star_R1 = FR1 + SR * (U_star_R1 - rhoR)
    F_star_R2 = FR2 + SR * (U_star_R2 - rhoR * uR)
    F_star_R3 = FR3 + SR * (U_star_R3 - rhoR * ER)

    # Flux assembly based on wave speeds
    f1 = np.where(SL >= 0.0, FL1,
         np.where(S_star >= 0.0, F_star_L1,
         np.where(SR >= 0.0, F_star_R1, FR1)))

    f2 = np.where(SL >= 0.0, FL2,
         np.where(S_star >= 0.0, F_star_L2,
         np.where(SR >= 0.0, F_star_R2, FR2)))

    f3 = np.where(SL >= 0.0, FL3,
         np.where(S_star >= 0.0, F_star_L3,
         np.where(SR >= 0.0, F_star_R3, FR3)))

    return np.column_stack((f1, f2, f3))
