"""
Numba JIT compiled computational kernels for pyshockflow.
Provides ultra-fast, single-pass fused loops for numerical fluxes, MUSCL reconstruction,
slope limiters, source terms, and time integration.
"""

import numpy as np

try:
    from numba import njit
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False
    # No-op decorator fallback if numba is not installed
    def njit(*args, **kwargs):
        def decorator(func):
            return func
        return decorator


# ---------------------------------------------------------------------------
# Slope Limiters (Numba inlined)
# ---------------------------------------------------------------------------

@njit(fastmath=True, inline='always')
def _limiter_van_albada(r):
    return (r * r + r) / (1.0 + r * r)

@njit(fastmath=True, inline='always')
def _limiter_van_leer(r):
    ar = abs(r)
    return (r + ar) / (1.0 + ar)

@njit(fastmath=True, inline='always')
def _limiter_minmod(r):
    return max(0.0, min(1.0, r))

@njit(fastmath=True, inline='always')
def _limiter_superbee(r):
    return max(0.0, max(min(2.0 * r, 1.0), min(r, 2.0)))

@njit(fastmath=True, inline='always')
def _evaluate_limiter(r, limiter_code):
    # limiter_code: 0=van albada, 1=van leer, 2=min-mod, 3=superbee, 4=none
    if limiter_code == 0:
        return _limiter_van_albada(r)
    elif limiter_code == 1:
        return _limiter_van_leer(r)
    elif limiter_code == 2:
        return _limiter_minmod(r)
    elif limiter_code == 3:
        return _limiter_superbee(r)
    else:
        return 1.0


# ---------------------------------------------------------------------------
# MUSCL Reconstruction (Numba)
# ---------------------------------------------------------------------------

@njit(fastmath=True)
def apply_muscl_reconstruction_numba(rho, u, p, x, ms, me,
                                     rhoL, uL, pL, rhoR, uR, pR,
                                     limiter_code):
    """Apply MUSCL high-order reconstruction to interior faces in-place."""
    for i in range(ms, me + 1):
        # Stencil points: il-1, il, ir, ir+1
        # In halo-inclusive indexing: il=i, ir=i+1
        # il-1 = i-1, il = i, ir = i+1, ir+1 = i+2
        
        dx_lm_l = x[i]     - x[i - 1]
        dx_l_r  = x[i + 1] - x[i]
        dx_r_rp = x[i + 2] - x[i + 1]
        
        # 1. Density
        r_l_rho  = ((rho[i] - rho[i - 1]) / dx_lm_l) / ((rho[i + 1] - rho[i]) / dx_l_r + 1e-6)
        r_r_rho  = ((rho[i + 1] - rho[i]) / dx_l_r)  / ((rho[i + 2] - rho[i + 1]) / dx_r_rp + 1e-6)
        psi_l_rho = _evaluate_limiter(r_l_rho, limiter_code)
        psi_r_rho = _evaluate_limiter(r_r_rho, limiter_code)
        rhoL[i] = rho[i]     + 0.5 * psi_l_rho * (rho[i + 1] - rho[i])
        rhoR[i] = rho[i + 1] - 0.5 * psi_r_rho * (rho[i + 2] - rho[i + 1])
        
        # 2. Velocity
        r_l_u  = ((u[i] - u[i - 1]) / dx_lm_l) / ((u[i + 1] - u[i]) / dx_l_r + 1e-6)
        r_r_u  = ((u[i + 1] - u[i]) / dx_l_r)  / ((u[i + 2] - u[i + 1]) / dx_r_rp + 1e-6)
        psi_l_u = _evaluate_limiter(r_l_u, limiter_code)
        psi_r_u = _evaluate_limiter(r_r_u, limiter_code)
        uL[i] = u[i]     + 0.5 * psi_l_u * (u[i + 1] - u[i])
        uR[i] = u[i + 1] - 0.5 * psi_r_u * (u[i + 2] - u[i + 1])
        
        # 3. Pressure
        r_l_p  = ((p[i] - p[i - 1]) / dx_lm_l) / ((p[i + 1] - p[i]) / dx_l_r + 1e-6)
        r_r_p  = ((p[i + 1] - p[i]) / dx_l_r)  / ((p[i + 2] - p[i + 1]) / dx_r_rp + 1e-6)
        psi_l_p = _evaluate_limiter(r_l_p, limiter_code)
        psi_r_p = _evaluate_limiter(r_r_p, limiter_code)
        pL[i] = p[i]     + 0.5 * psi_l_p * (p[i + 1] - p[i])
        pR[i] = p[i + 1] - 0.5 * psi_r_p * (p[i + 2] - p[i + 1])


# ---------------------------------------------------------------------------
# Roe Numerical Flux: Ideal Gas (Toro formulation)
# ---------------------------------------------------------------------------

@njit(fastmath=True)
def compute_roe_flux_ideal_numba(rhoL, rhoR, uL, uR, pL, pR, gmma,
                                 entropy_fix_active, fix_coeff):
    n = len(rhoL)
    flux = np.empty((n, 3))
    gm1 = gmma - 1.0
    
    for i in range(n):
        rL = rhoL[i]
        rR = rhoR[i]
        vL = uL[i]
        vR = uR[i]
        prL = pL[i]
        prR = pR[i]
        
        eL = prL / (gm1 * rL)
        eR = prR / (gm1 * rR)
        
        htL = 0.5 * vL * vL + eL + prL / rL
        htR = 0.5 * vR * vR + eR + prR / rR
        
        sqrL = np.sqrt(rL)
        sqrR = np.sqrt(rR)
        denom = sqrL + sqrR
        
        uAVG = (sqrL * vL + sqrR * vR) / denom
        hAVG = (sqrL * htL + sqrR * htR) / denom
        aAVG = np.sqrt(max(1e-10, gm1 * (hAVG - 0.5 * uAVG * uAVG)))
        rhoAVG = sqrL * sqrR
        
        dp = prR - prL
        du = vR - vL
        drho = rR - rL
        a2inv = 1.0 / (aAVG * aAVG)
        
        alpha1 = 0.5 * a2inv * (dp - rhoAVG * aAVG * du)
        alpha2 = drho - dp * a2inv
        alpha3 = 0.5 * a2inv * (dp + rhoAVG * aAVG * du)
        
        lam1 = uAVG - aAVG
        lam2 = uAVG
        lam3 = uAVG + aAVG
        
        alam1 = abs(lam1)
        alam2 = abs(lam2)
        alam3 = abs(lam3)
        
        if entropy_fix_active and fix_coeff > 0.0:
            delta = fix_coeff * aAVG
            if alam1 < delta:
                alam1 = 0.5 * (alam1 * alam1 / delta + delta)
            if alam2 < delta:
                alam2 = 0.5 * (alam2 * alam2 / delta + delta)
            if alam3 < delta:
                alam3 = 0.5 * (alam3 * alam3 / delta + delta)
                
        fl1 = rL * vL
        fl2 = rL * vL * vL + prL
        fl3 = vL * (rL * (0.5 * vL * vL + eL) + prL)
        
        fr1 = rR * vR
        fr2 = rR * vR * vR + prR
        fr3 = vR * (rR * (0.5 * vR * vR + eR) + prR)
        
        diss1 = (alam1 * alpha1 + 
                 alam2 * alpha2 + 
                 alam3 * alpha3)
                 
        diss2 = (alam1 * alpha1 * (uAVG - aAVG) + 
                 alam2 * alpha2 * uAVG + 
                 alam3 * alpha3 * (uAVG + aAVG))
                 
        diss3 = (alam1 * alpha1 * (hAVG - uAVG * aAVG) + 
                 alam2 * alpha2 * (0.5 * uAVG * uAVG) + 
                 alam3 * alpha3 * (hAVG + uAVG * aAVG))
                 
        flux[i, 0] = 0.5 * (fl1 + fr1 - diss1)
        flux[i, 1] = 0.5 * (fl2 + fr2 - diss2)
        flux[i, 2] = 0.5 * (fl3 + fr3 - diss3)
        
    return flux


# ---------------------------------------------------------------------------
# Roe Numerical Flux: Real Gas (Arabi et al. 2017)
# ---------------------------------------------------------------------------

@njit(fastmath=True)
def compute_roe_flux_arabi_numba(rhoL, rhoR, uL, uR, pL, pR,
                                 eL, eR, aL, aR,
                                 entropy_fix_active, fix_coeff):
    n = len(rhoL)
    flux = np.empty((n, 3))
    
    for i in range(n):
        rL = rhoL[i]
        rR = rhoR[i]
        vL = uL[i]
        vR = uR[i]
        prL = pL[i]
        prR = pR[i]
        enL = eL[i]
        enR = eR[i]
        asL = aL[i]
        asR = aR[i]
        
        htL = 0.5 * vL * vL + enL + prL / rL
        htR = 0.5 * vR * vR + enR + prR / rR
        
        sqrL = np.sqrt(rL)
        sqrR = np.sqrt(rR)
        denom = sqrL + sqrR
        
        uAVG = (sqrL * vL + sqrR * vR) / denom
        hAVG = (sqrL * htL + sqrR * htR) / denom
        aAVG = (sqrL * asL + sqrR * asR) / denom
        rhoAVG = sqrL * sqrR
        
        dp = prR - prL
        du = vR - vL
        drho = rR - rL
        a2inv = 1.0 / (aAVG * aAVG)
        
        # Arabi eigenvalues: u+a, u-a, u
        lam1 = uAVG + aAVG
        lam2 = uAVG - aAVG
        lam3 = uAVG
        
        alam1 = abs(lam1)
        alam2 = abs(lam2)
        alam3 = abs(lam3)
        
        if entropy_fix_active and fix_coeff > 0.0:
            delta = fix_coeff * aAVG
            if alam1 < delta:
                alam1 = 0.5 * (alam1 * alam1 / delta + delta)
            if alam2 < delta:
                alam2 = 0.5 * (alam2 * alam2 / delta + delta)
            if alam3 < delta:
                alam3 = 0.5 * (alam3 * alam3 / delta + delta)
                
        # Arabi wave strengths
        alpha1 = 0.5 * a2inv * (dp + rhoAVG * aAVG * du)
        alpha2 = 0.5 * a2inv * (dp - rhoAVG * aAVG * du)
        alpha3 = drho - dp * a2inv
        
        # Physical Euler fluxes
        fL1 = rL * vL
        fL2 = rL * vL * vL + prL
        fL3 = vL * (rL * (0.5 * vL * vL + enL) + prL)
        
        fr1 = rR * vR
        fr2 = rR * vR * vR + prR
        fr3 = vR * (rR * (0.5 * vR * vR + enR) + prR)
        
        # Dissipation components
        d1 = alam1 * alpha1 + alam2 * alpha2 + alam3 * alpha3
        d2 = (lam1 * alam1 * alpha1 +
              lam2 * alam2 * alpha2 +
              uAVG * alam3 * alpha3)
              
        # Arabi X energy term
        X = ((rR * vR * htR - rL * vL * htL) -
             (hAVG + uAVG * aAVG) * lam1 * alpha1 -
             (hAVG - uAVG * aAVG) * lam2 * alpha2)
        if uAVG < 0:
            X = -X
            
        d3 = ((hAVG + uAVG * aAVG) * alam1 * alpha1 +
              (hAVG - uAVG * aAVG) * alam2 * alpha2 + X)
              
        flux[i, 0] = 0.5 * (fL1 + fr1 - d1)
        flux[i, 1] = 0.5 * (fL2 + fr2 - d2)
        flux[i, 2] = 0.5 * (fL3 + fr3 - d3)
        
    return flux


# ---------------------------------------------------------------------------
# Roe Numerical Flux: Real Gas (Vinokur & Montagné 1990)
# ---------------------------------------------------------------------------

@njit(fastmath=True)
def compute_roe_flux_vinokur_numba(rhoL, rhoR, uL, uR, pL, pR,
                                   eL, eR,
                                   chiL, chiR, chiM,
                                   kappaL, kappaR, kappaM,
                                   entropy_fix_active, fix_coeff):
    n = len(rhoL)
    flux = np.empty((n, 3))
    
    for i in range(n):
        rL = rhoL[i]
        rR = rhoR[i]
        vL = uL[i]
        vR = uR[i]
        prL = pL[i]
        prR = pR[i]
        enL = eL[i]
        enR = eR[i]
        
        htL = 0.5 * vL * vL + enL + prL / rL
        htR = 0.5 * vR * vR + enR + prR / rR
        
        sqrL = np.sqrt(rL)
        sqrR = np.sqrt(rR)
        alpha = sqrL / (sqrL + sqrR)
        
        du = vR - vL
        dp = prR - prL
        drho = rR - rL
        
        uAVG = alpha * vL + (1.0 - alpha) * vR
        htAVG = alpha * htL + (1.0 - alpha) * htR
        hL = htL - 0.5 * vL * vL
        hR = htR - 0.5 * vR * vR
        hAVG = (alpha * hL + (1.0 - alpha) * hR +
                0.5 * alpha * (1.0 - alpha) * du * du)
                
        rhoeL = rL * enL
        rhoeR = rR * enR
        delta_rhoe = rhoeR - rhoeL
        
        chiHat = (chiL[i] + chiR[i] + 4.0 * chiM[i]) / 6.0
        kappaHat = (kappaL[i] + kappaR[i] + 4.0 * kappaM[i]) / 6.0
        
        error_term = dp - chiHat * drho - kappaHat * delta_rhoe
        hM = 0.5 * (hL + hR)
        csquare_L = chiL[i] + kappaL[i] * hL
        csquare_R = chiR[i] + kappaR[i] * hR
        csquare_M = chiM[i] + kappaM[i] * hM
        sHat = (csquare_L + csquare_R + 4.0 * csquare_M) / 6.0
        D_term = (sHat * drho) * (sHat * drho) + dp * dp
        
        denom_proj = D_term - dp * error_term
        if denom_proj == 0.0:
            safe_denom = 1.0
        else:
            safe_denom = denom_proj
            
        if drho == 0.0:
            chiAVG = chiHat
        else:
            chiAVG = (D_term * chiHat + sHat * sHat * drho * error_term) / safe_denom
            
        if dp == 0.0:
            kappaAVG = kappaHat
        else:
            kappaAVG = (D_term * kappaHat) / safe_denom
            
        aAVG = np.sqrt(max(1e-10, chiAVG + kappaAVG * hAVG))
        
        # Conservative variable jumps
        u3L = rL * (0.5 * vL * vL + enL)
        u3R = rR * (0.5 * vR * vR + enR)
        dU0 = drho
        dU1 = rR * vR - rL * vL
        dU2 = u3R - u3L
        
        k1 = 0.5 * kappaAVG * uAVG * uAVG + kappaAVG
        k2 = 0.5 * uAVG * uAVG - chiAVG / max(1e-10, kappaAVG)
        
        # Vinokur eigenvalues: u, u+a, u-a
        lam1 = uAVG
        lam2 = uAVG + aAVG
        lam3 = uAVG - aAVG
        
        alam1 = abs(lam1)
        alam2 = abs(lam2)
        alam3 = abs(lam3)
        
        if entropy_fix_active and fix_coeff > 0.0:
            delta = fix_coeff * aAVG
            if alam1 < delta:
                alam1 = 0.5 * (alam1 * alam1 / delta + delta)
            if alam2 < delta:
                alam2 = 0.5 * (alam2 * alam2 / delta + delta)
            if alam3 < delta:
                alam3 = 0.5 * (alam3 * alam3 / delta + delta)
                
        a2inv = 1.0 / (aAVG * aAVG)
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
              
        aw0 = alam1 * w0
        aw1 = alam2 * w1
        aw2 = alam3 * w2
        
        d1 = aw0 + aw1 + aw2
        d2 = uAVG * aw0 + (uAVG + aAVG) * aw1 + (uAVG - aAVG) * aw2
        d3 = k2 * aw0 + (htAVG + aAVG * uAVG) * aw1 + (htAVG - aAVG * uAVG) * aw2
        
        fL1 = rL * vL
        fL2 = rL * vL * vL + prL
        fL3 = vL * (u3L + prL)
        
        fr1 = rR * vR
        fr2 = rR * vR * vR + prR
        fr3 = vR * (u3R + prR)
        
        flux[i, 0] = 0.5 * (fL1 + fr1 - d1)
        flux[i, 1] = 0.5 * (fL2 + fr2 - d2)
        flux[i, 2] = 0.5 * (fL3 + fr3 - d3)
        
    return flux


# ---------------------------------------------------------------------------
# HLLC Numerical Flux: Ideal Gas (Numba)
# ---------------------------------------------------------------------------

@njit(fastmath=True)
def compute_hllc_flux_ideal_numba(rhoL, rhoR, uL, uR, pL, pR, gmma):
    n = len(rhoL)
    flux = np.empty((n, 3))
    gm1 = gmma - 1.0
    
    for i in range(n):
        rL = rhoL[i]
        rR = rhoR[i]
        vL = uL[i]
        vR = uR[i]
        prL = pL[i]
        prR = pR[i]
        
        enL = prL / (gm1 * rL)
        enR = prR / (gm1 * rR)
        aL = np.sqrt(gmma * prL / rL)
        aR = np.sqrt(gmma * prR / rR)
        
        EL = enL + 0.5 * vL * vL
        ER = enR + 0.5 * vR * vR
        
        fL1 = rL * vL
        fL2 = rL * vL * vL + prL
        fL3 = vL * (rL * EL + prL)
        
        fr1 = rR * vR
        fr2 = rR * vR * vR + prR
        fr3 = vR * (rR * ER + prR)
        
        SL = min(vL - aL, vR - aR)
        SR = max(vL + aL, vR + aR)
        
        num = prR - prL + rL * vL * (SL - vL) - rR * vR * (SR - vR)
        den = rL * (SL - vL) - rR * (SR - vR)
        
        if abs(den) < 1e-14:
            S_star = 0.5 * (vL + vR)
        else:
            S_star = num / den
            
        if SL >= 0.0:
            flux[i, 0] = fL1
            flux[i, 1] = fL2
            flux[i, 2] = fL3
        elif S_star >= 0.0:
            diff_L = SL - S_star
            safe_diff_L = 1e-30 if abs(diff_L) < 1e-30 else diff_L
            diff_SL_uL = SL - vL
            safe_SL_uL = 1e-30 if abs(diff_SL_uL) < 1e-30 else diff_SL_uL
            
            factor_L = rL * diff_SL_uL / safe_diff_L
            E_star_L = EL + (S_star - vL) * (S_star + prL / (rL * safe_SL_uL))
            
            u_star_L1 = factor_L
            u_star_L2 = factor_L * S_star
            u_star_L3 = factor_L * E_star_L
            
            flux[i, 0] = fL1 + SL * (u_star_L1 - rL)
            flux[i, 1] = fL2 + SL * (u_star_L2 - rL * vL)
            flux[i, 2] = fL3 + SL * (u_star_L3 - rL * EL)
        elif SR >= 0.0:
            diff_R = SR - S_star
            safe_diff_R = 1e-30 if abs(diff_R) < 1e-30 else diff_R
            diff_SR_uR = SR - vR
            safe_SR_uR = 1e-30 if abs(diff_SR_uR) < 1e-30 else diff_SR_uR
            
            factor_R = rR * diff_SR_uR / safe_diff_R
            E_star_R = ER + (S_star - vR) * (S_star + prR / (rR * safe_SR_uR))
            
            u_star_R1 = factor_R
            u_star_R2 = factor_R * S_star
            u_star_R3 = factor_R * E_star_R
            
            flux[i, 0] = fr1 + SR * (u_star_R1 - rR)
            flux[i, 1] = fr2 + SR * (u_star_R2 - rR * vR)
            flux[i, 2] = fr3 + SR * (u_star_R3 - rR * ER)
        else:
            flux[i, 0] = fr1
            flux[i, 1] = fr2
            flux[i, 2] = fr3
            
    return flux


# ---------------------------------------------------------------------------
# HLLC Numerical Flux: Real Gas (Numba)
# ---------------------------------------------------------------------------

@njit(fastmath=True)
def compute_hllc_flux_real_numba(rhoL, rhoR, uL, uR, pL, pR,
                                 eL, eR, aL_arr, aR_arr):
    n = len(rhoL)
    flux = np.empty((n, 3))
    
    for i in range(n):
        rL = rhoL[i]
        rR = rhoR[i]
        vL = uL[i]
        vR = uR[i]
        prL = pL[i]
        prR = pR[i]
        enL = eL[i]
        enR = eR[i]
        aL = aL_arr[i]
        aR = aR_arr[i]
        
        EL = enL + 0.5 * vL * vL
        ER = enR + 0.5 * vR * vR
        
        fL1 = rL * vL
        fL2 = rL * vL * vL + prL
        fL3 = vL * (rL * EL + prL)
        
        fr1 = rR * vR
        fr2 = rR * vR * vR + prR
        fr3 = vR * (rR * ER + prR)
        
        SL = min(vL - aL, vR - aR)
        SR = max(vL + aL, vR + aR)
        
        num = prR - prL + rL * vL * (SL - vL) - rR * vR * (SR - vR)
        den = rL * (SL - vL) - rR * (SR - vR)
        
        if abs(den) < 1e-14:
            S_star = 0.5 * (vL + vR)
        else:
            S_star = num / den
            
        if SL >= 0.0:
            flux[i, 0] = fL1
            flux[i, 1] = fL2
            flux[i, 2] = fL3
        elif S_star >= 0.0:
            diff_L = SL - S_star
            safe_diff_L = 1e-30 if abs(diff_L) < 1e-30 else diff_L
            diff_SL_uL = SL - vL
            safe_SL_uL = 1e-30 if abs(diff_SL_uL) < 1e-30 else diff_SL_uL
            
            factor_L = rL * diff_SL_uL / safe_diff_L
            E_star_L = EL + (S_star - vL) * (S_star + prL / (rL * safe_SL_uL))
            
            u_star_L1 = factor_L
            u_star_L2 = factor_L * S_star
            u_star_L3 = factor_L * E_star_L
            
            flux[i, 0] = fL1 + SL * (u_star_L1 - rL)
            flux[i, 1] = fL2 + SL * (u_star_L2 - rL * vL)
            flux[i, 2] = fL3 + SL * (u_star_L3 - rL * EL)
        elif SR >= 0.0:
            diff_R = SR - S_star
            safe_diff_R = 1e-30 if abs(diff_R) < 1e-30 else diff_R
            diff_SR_uR = SR - vR
            safe_SR_uR = 1e-30 if abs(diff_SR_uR) < 1e-30 else diff_SR_uR
            
            factor_R = rR * diff_SR_uR / safe_diff_R
            E_star_R = ER + (S_star - vR) * (S_star + prR / (rR * safe_SR_uR))
            
            u_star_R1 = factor_R
            u_star_R2 = factor_R * S_star
            u_star_R3 = factor_R * E_star_R
            
            flux[i, 0] = fr1 + SR * (u_star_R1 - rR)
            flux[i, 1] = fr2 + SR * (u_star_R2 - rR * vR)
            flux[i, 2] = fr3 + SR * (u_star_R3 - rR * ER)
        else:
            flux[i, 0] = fr1
            flux[i, 1] = fr2
            flux[i, 2] = fr3
            
    return flux


# ---------------------------------------------------------------------------
# Residual Assembly & Solution Update (Numba)
# ---------------------------------------------------------------------------

@njit(fastmath=True)
def compute_residuals_numba(flux, source, dV, dt, areaReference):
    """Fused residual assembly across all physical interior cells."""
    n = len(dV) - 2
    residuals = np.empty((n, 3))
    
    for i in range(n):
        inv_dv = dt / dV[i + 1]
        
        # Flux divergence
        df1 = (flux[i, 0] - flux[i + 1, 0]) * areaReference
        df2 = (flux[i, 1] - flux[i + 1, 1]) * areaReference
        df3 = (flux[i, 2] - flux[i + 1, 2]) * areaReference
        
        # Add volumetric source terms
        s1 = source[i + 1, 0] * dV[i + 1]
        s2 = source[i + 1, 1] * dV[i + 1]
        s3 = source[i + 1, 2] * dV[i + 1]
        
        residuals[i, 0] = inv_dv * (df1 + s1)
        residuals[i, 1] = inv_dv * (df2 + s2)
        residuals[i, 2] = inv_dv * (df3 + s3)
        
    return residuals


@njit(fastmath=True)
def update_solution_numba(u1, u2, u3, residuals):
    """In-place update of conservative variables."""
    n = len(residuals)
    for i in range(n):
        u1[i + 1] += residuals[i, 0]
        u2[i + 1] += residuals[i, 1]
        u3[i + 1] += residuals[i, 2]


@njit(fastmath=True)
def primitives_from_conservatives_ideal_numba(u1, u2, u3, gmma, rho, u, p, e):
    """Recover ideal gas primitive variables directly in compiled machine code."""
    n = len(u1) - 2
    gm1 = gmma - 1.0
    for i in range(1, n + 1):
        r = u1[i]
        v = u2[i] / r
        en = (u3[i] / r) - 0.5 * v * v
        pr = gm1 * r * en
        
        rho[i] = r
        u[i] = v
        p[i] = pr
        e[i] = en
