"""
Vectorized Roe approximate Riemann solvers for the 1D Euler equations.

These functions compute numerical fluxes at ALL interfaces simultaneously
using NumPy array operations, replacing the per-interface class instantiation
and Python loops in the original scalar implementation.

Three variants are provided:
- ``compute_roe_flux_ideal``: Standard Roe (Toro) for ideal gas
- ``compute_roe_flux_arabi``: Arabi et al. (2017) for real gases
- ``compute_roe_flux_vinokur``: Vinokur & Montagné (1990) for real gases
"""

import numpy as np


# ---------------------------------------------------------------------------
# Vectorized Harten entropy fix
# ---------------------------------------------------------------------------

def apply_entropy_fix_vec(lam, a_avg, kappa):
    """Harten entropy fix applied element-wise to an eigenvalue array.

    Parameters
    ----------
    lam : ndarray (N,)
        Raw eigenvalue at each interface.
    a_avg : ndarray (N,)
        Roe-averaged sound speed at each interface.
    kappa : float
        Fix coefficient (typically 0.2).

    Returns
    -------
    ndarray (N,)
        Corrected absolute eigenvalues.
    """
    delta = kappa * a_avg
    abs_lam = np.abs(lam)
    return np.where(abs_lam < delta,
                    0.5 * (lam ** 2 / delta + delta),
                    abs_lam)


def _abs_eigenvalues(lam1, lam2, lam3, a_avg, entropy_fix_active, fix_coeff):
    """Return |eigenvalues| with optional entropy fix."""
    if entropy_fix_active:
        return (apply_entropy_fix_vec(lam1, a_avg, fix_coeff),
                apply_entropy_fix_vec(lam2, a_avg, fix_coeff),
                apply_entropy_fix_vec(lam3, a_avg, fix_coeff))
    return np.abs(lam1), np.abs(lam2), np.abs(lam3)


# ---------------------------------------------------------------------------
# Physical Euler flux from primitives
# ---------------------------------------------------------------------------

def _euler_flux(rho, u, p, e):
    """Euler flux vector from primitive variables.

    Returns f1, f2, f3 each of shape (N,).
    """
    et = e + 0.5 * u ** 2
    f1 = rho * u
    f2 = rho * u ** 2 + p
    f3 = u * (rho * et + p)
    return f1, f2, f3


# ---------------------------------------------------------------------------
# Standard Roe flux — ideal gas (Toro formulation)
# ---------------------------------------------------------------------------

def compute_roe_flux_ideal(rhoL, rhoR, uL, uR, pL, pR, gmma,
                           entropy_fix_active, fix_coeff):
    """Vectorized Roe flux for ideal gas.

    Parameters
    ----------
    rhoL, rhoR, uL, uR, pL, pR : ndarray (N,)
        Left / right primitive states at each of the *N* interfaces.
    gmma : float
        Ratio of specific heats γ.
    entropy_fix_active : bool
    fix_coeff : float

    Returns
    -------
    flux : ndarray (N, 3)
        Numerical flux at each interface.
    """
    gm1 = gmma - 1.0

    # Static energy (ideal gas)
    eL = pL / (gm1 * rhoL)
    eR = pR / (gm1 * rhoR)

    # Total enthalpy
    htL = 0.5 * uL ** 2 + eL + pL / rhoL
    htR = 0.5 * uR ** 2 + eR + pR / rhoR

    # Roe averages
    sqrL = np.sqrt(rhoL)
    sqrR = np.sqrt(rhoR)
    denom = sqrL + sqrR
    uAVG = (sqrL * uL + sqrR * uR) / denom
    hAVG = (sqrL * htL + sqrR * htR) / denom
    aAVG = np.sqrt(gm1 * (hAVG - 0.5 * uAVG ** 2))
    rhoAVG = sqrL * sqrR

    # Jumps
    dp = pR - pL
    du = uR - uL
    drho = rhoR - rhoL
    a2inv = 1.0 / (aAVG ** 2)

    # Wave strengths (Toro ordering: u-a, u, u+a)
    alpha1 = 0.5 * a2inv * (dp - rhoAVG * aAVG * du)
    alpha2 = drho - dp * a2inv
    alpha3 = 0.5 * a2inv * (dp + rhoAVG * aAVG * du)

    # Eigenvalues
    lam1 = uAVG - aAVG
    lam2 = uAVG
    lam3 = uAVG + aAVG

    absLam1, absLam2, absLam3 = _abs_eigenvalues(
        lam1, lam2, lam3, aAVG, entropy_fix_active, fix_coeff)

    # Physical Euler fluxes
    fL1, fL2, fL3 = _euler_flux(rhoL, uL, pL, eL)
    fR1, fR2, fR3 = _euler_flux(rhoR, uR, pR, eR)

    # Dissipation = eigenvector_matrix @ diag(|λ|) @ wave_strengths
    #   expanded algebraically to avoid (N,3,3) allocations
    d1 = (alpha1 * absLam1 +
          alpha2 * absLam2 +
          alpha3 * absLam3)
    d2 = (alpha1 * absLam1 * (uAVG - aAVG) +
          alpha2 * absLam2 * uAVG +
          alpha3 * absLam3 * (uAVG + aAVG))
    d3 = (alpha1 * absLam1 * (hAVG - uAVG * aAVG) +
          alpha2 * absLam2 * 0.5 * uAVG ** 2 +
          alpha3 * absLam3 * (hAVG + uAVG * aAVG))

    # Assemble Roe flux
    flux = np.empty((len(rhoL), 3))
    flux[:, 0] = 0.5 * (fL1 + fR1) - 0.5 * d1
    flux[:, 1] = 0.5 * (fL2 + fR2) - 0.5 * d2
    flux[:, 2] = 0.5 * (fL3 + fR3) - 0.5 * d3
    return flux


# ---------------------------------------------------------------------------
# Arabi real-gas Roe flux
# ---------------------------------------------------------------------------

def compute_roe_flux_arabi(rhoL, rhoR, uL, uR, pL, pR,
                           eL, eR, aL, aR,
                           entropy_fix_active, fix_coeff):
    """Vectorized Arabi et al. (2017) real-gas Roe flux.

    Parameters
    ----------
    rhoL, rhoR, uL, uR, pL, pR : ndarray (N,)
    eL, eR : ndarray (N,)
        Static (internal) energy.
    aL, aR : ndarray (N,)
        Speed of sound.
    entropy_fix_active : bool
    fix_coeff : float

    Returns
    -------
    flux : ndarray (N, 3)
    """
    # Total enthalpy
    htL = 0.5 * uL ** 2 + eL + pL / rhoL
    htR = 0.5 * uR ** 2 + eR + pR / rhoR

    # Roe averages
    sqrL = np.sqrt(rhoL)
    sqrR = np.sqrt(rhoR)
    denom = sqrL + sqrR
    uAVG = (sqrL * uL + sqrR * uR) / denom
    hAVG = (sqrL * htL + sqrR * htR) / denom
    aAVG = (sqrL * aL + sqrR * aR) / denom   # Arabi: directly averaged
    rhoAVG = sqrL * sqrR

    # Jumps
    dp = pR - pL
    du = uR - uL
    drho = rhoR - rhoL
    a2inv = 1.0 / (aAVG ** 2)

    # Eigenvalues (Arabi ordering: u+a, u-a, u)
    lam1 = uAVG + aAVG
    lam2 = uAVG - aAVG
    lam3 = uAVG

    absLam1, absLam2, absLam3 = _abs_eigenvalues(
        lam1, lam2, lam3, aAVG, entropy_fix_active, fix_coeff)

    # Wave strengths (Arabi formulation)
    alpha1 = 0.5 * a2inv * (dp + rhoAVG * aAVG * du)
    alpha2 = 0.5 * a2inv * (dp - rhoAVG * aAVG * du)
    alpha3 = drho - dp * a2inv

    # Physical Euler fluxes
    fL1, fL2, fL3 = _euler_flux(rhoL, uL, pL, eL)
    fR1, fR2, fR3 = _euler_flux(rhoR, uR, pR, eR)

    # Dissipation components
    d1 = absLam1 * alpha1 + absLam2 * alpha2 + absLam3 * alpha3
    d2 = (lam1 * absLam1 * alpha1 +
          lam2 * absLam2 * alpha2 +
          uAVG * absLam3 * alpha3)

    # X term for energy equation (Arabi-specific)
    X = ((rhoR * uR * htR - rhoL * uL * htL) -
         (hAVG + uAVG * aAVG) * lam1 * alpha1 -
         (hAVG - uAVG * aAVG) * lam2 * alpha2)
    X = np.where(uAVG >= 0, X, -X)

    d3 = ((hAVG + uAVG * aAVG) * absLam1 * alpha1 +
          (hAVG - uAVG * aAVG) * absLam2 * alpha2 + X)

    # Assemble
    flux = np.empty((len(rhoL), 3))
    flux[:, 0] = 0.5 * (fL1 + fR1) - 0.5 * d1
    flux[:, 1] = 0.5 * (fL2 + fR2) - 0.5 * d2
    flux[:, 2] = 0.5 * (fL3 + fR3) - 0.5 * d3
    return flux


# ---------------------------------------------------------------------------
# Vinokur–Montagné real-gas Roe flux
# ---------------------------------------------------------------------------

def compute_roe_flux_vinokur(rhoL, rhoR, uL, uR, pL, pR,
                             eL, eR,
                             chiL, chiR, chiM,
                             kappaL, kappaR, kappaM,
                             entropy_fix_active, fix_coeff):
    """Vectorized Vinokur & Montagné (1990) real-gas Roe flux.

    Parameters
    ----------
    rhoL, rhoR, uL, uR, pL, pR : ndarray (N,)
    eL, eR : ndarray (N,)
        Static (internal) energy.
    chiL, chiR, chiM : ndarray (N,)
        Pressure derivative χ = (∂p/∂ρ)|e − (e/ρ)(∂p/∂e)|ρ for left, right,
        and arithmetic-mean states.
    kappaL, kappaR, kappaM : ndarray (N,)
        Pressure derivative κ = (1/ρ)(∂p/∂e)|ρ for left, right, and mean.
    entropy_fix_active : bool
    fix_coeff : float

    Returns
    -------
    flux : ndarray (N, 3)
    """
    # Total enthalpy
    htL = 0.5 * uL ** 2 + eL + pL / rhoL
    htR = 0.5 * uR ** 2 + eR + pR / rhoR

    # Roe averaging weight
    sqrL = np.sqrt(rhoL)
    sqrR = np.sqrt(rhoR)
    alpha = sqrL / (sqrL + sqrR)

    # Jumps
    du = uR - uL
    dp = pR - pL
    drho = rhoR - rhoL

    # Averaged variables (Vinokur formulation)
    uAVG = alpha * uL + (1.0 - alpha) * uR
    htAVG = alpha * htL + (1.0 - alpha) * htR
    hL = htL - 0.5 * uL ** 2
    hR = htR - 0.5 * uR ** 2
    hAVG = (alpha * hL + (1.0 - alpha) * hR +
            0.5 * alpha * (1.0 - alpha) * du ** 2)

    # Mean-state quantities for projection
    rhoeL = rhoL * eL
    rhoeR = rhoR * eR
    delta_rhoe = rhoeR - rhoeL

    # Simpson-like averaging of chi, kappa
    chiHat = (chiL + chiR + 4.0 * chiM) / 6.0
    kappaHat = (kappaL + kappaR + 4.0 * kappaM) / 6.0

    # Projection procedure
    error_term = dp - chiHat * drho - kappaHat * delta_rhoe
    hM = 0.5 * (hL + hR)
    csquare_L = chiL + kappaL * hL
    csquare_R = chiR + kappaR * hR
    csquare_M = chiM + kappaM * hM
    sHat = (csquare_L + csquare_R + 4.0 * csquare_M) / 6.0
    D_term = (sHat * drho) ** 2 + dp ** 2

    denom_proj = D_term - dp * error_term
    safe_denom = np.where(denom_proj == 0.0, 1.0, denom_proj)
    chiAVG = np.where(
        drho == 0,
        chiHat,
        (D_term * chiHat + sHat ** 2 * drho * error_term) / safe_denom)
    kappaAVG = np.where(
        dp == 0,
        kappaHat,
        (D_term * kappaHat) / safe_denom)

    aAVG = np.sqrt(chiAVG + kappaAVG * hAVG)

    # Conservative-variable jumps
    u3L = rhoL * (0.5 * uL ** 2 + eL)
    u3R = rhoR * (0.5 * uR ** 2 + eR)
    dU0 = drho                      # Δ(ρ)
    dU1 = rhoR * uR - rhoL * uL    # Δ(ρu)
    dU2 = u3R - u3L                 # Δ(ρe_t)

    # Auxiliary coefficients (from eigenvector matrices)
    k1 = 0.5 * kappaAVG * uAVG ** 2 + kappaAVG
    k2 = 0.5 * uAVG ** 2 - chiAVG / kappaAVG

    # Eigenvalues (Vinokur ordering: u, u+a, u-a)
    lam1 = uAVG
    lam2 = uAVG + aAVG
    lam3 = uAVG - aAVG

    absLam1, absLam2, absLam3 = _abs_eigenvalues(
        lam1, lam2, lam3, aAVG, entropy_fix_active, fix_coeff)

    # Wave strengths w = R⁻¹ @ ΔU  (expanded algebraically)
    a2inv = 1.0 / (aAVG ** 2)
    ainv = 1.0 / aAVG
    ku_a2 = kappaAVG * uAVG * a2inv
    k_a2 = kappaAVG * a2inv
    k1_a2 = k1 * a2inv
    u_ainv = uAVG * ainv

    w0 = (1.0 - k1_a2) * dU0 + ku_a2 * dU1 - k_a2 * dU2
    w1 = (0.5 * (k1_a2 - u_ainv) * dU0 -
          0.5 * (ku_a2 - ainv) * dU1 +
          0.5 * k_a2 * dU2)
    w2 = (0.5 * (k1_a2 + u_ainv) * dU0 -
          0.5 * (ku_a2 + ainv) * dU1 +
          0.5 * k_a2 * dU2)

    # |λ| * w
    aw0 = absLam1 * w0
    aw1 = absLam2 * w1
    aw2 = absLam3 * w2

    # Dissipation = R @ (|λ| * w), expanded
    d1 = aw0 + aw1 + aw2
    d2 = uAVG * aw0 + (uAVG + aAVG) * aw1 + (uAVG - aAVG) * aw2
    d3 = k2 * aw0 + (htAVG + aAVG * uAVG) * aw1 + (htAVG - aAVG * uAVG) * aw2

    # Physical Euler fluxes
    fL1, fL2, fL3 = _euler_flux(rhoL, uL, pL, eL)
    fR1, fR2, fR3 = _euler_flux(rhoR, uR, pR, eR)

    # Assemble
    flux = np.empty((len(rhoL), 3))
    flux[:, 0] = 0.5 * (fL1 + fR1) - 0.5 * d1
    flux[:, 1] = 0.5 * (fL2 + fR2) - 0.5 * d2
    flux[:, 2] = 0.5 * (fL3 + fR3) - 0.5 * d3
    return flux
