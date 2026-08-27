#!/usr/bin/env python3
"""
1-D Rayleigh-flow solver for a constant-area circular duct.

Assumptions
-----------
- Steady 1-D compressible flow
- Constant-area circular duct
- Perfect gas
- Constant gamma and R
- No wall friction
- No shaft work
- Heat transfer only
- Constant wall heat flux qw [W/m^2]

Boundary conditions
-------------------
Inlet:
    total pressure    pt_in
    total temperature Tt_in

Outlet:
    static pressure   p_out

For each prescribed wall heat flux qw, the inlet Mach number is
determined such that the integrated Rayleigh-flow solution satisfies
the prescribed outlet static pressure.

Author: ...
"""

import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import solve_ivp
from scipy.optimize import brentq


# =============================================================================
# USER INPUT
# =============================================================================

# Geometry
L = 10.0                  # duct length [m]
A = 0.01                  # cross-sectional area [m^2]

# Inlet stagnation conditions
pt_in = 101.325e3           # inlet total pressure [Pa]
Tt_in = 288.15             # inlet total temperature [K]

# Outlet static pressure
p_out = 90.0e3             # outlet static pressure [Pa]

# Wall heat fluxes
qw_values = np.array([
    0.0,
    1.0,
    5.0,
    10.0
]) * 1.0e3                # [W/m^2]

# Perfect-gas properties
gamma = 1.4
R = 287.06              # [J/(kg K)]

cp = gamma * R / (gamma - 1.0)


# =============================================================================
# NUMERICAL PARAMETERS
# =============================================================================

RTOL = 1.0e-5
ATOL = 1.0e-5
MAX_STEP = 0.02            # [m]

# Subsonic branch
M_MIN = 1.0e-5
M_MAX = 0.9999


# =============================================================================
# GEOMETRY
# =============================================================================

D = np.sqrt(4.0 * A / np.pi)
P = np.pi * D

print("=" * 75)
print("RAYLEIGH FLOW SOLVER")
print("=" * 75)

print(f"Duct length       : {L:.4f} m")
print(f"Cross-sectional A : {A:.6f} m^2")
print(f"Diameter          : {D:.6f} m")
print(f"Wetted perimeter  : {P:.6f} m")
print(f"Inlet total p     : {pt_in / 1e3:.3f} kPa")
print(f"Inlet total T     : {Tt_in:.3f} K")
print(f"Outlet static p   : {p_out / 1e3:.3f} kPa")
print("=" * 75)


# =============================================================================
# BASIC THERMODYNAMIC RELATIONS
# =============================================================================

def static_from_total(M, pt, Tt):
    """
    Convert stagnation properties to static properties
    for a perfect gas under isentropic stagnation relations.

    Parameters
    ----------
    M : float
        Mach number
    pt : float
        Total pressure [Pa]
    Tt : float
        Total temperature [K]

    Returns
    -------
    p : float
        Static pressure [Pa]
    T : float
        Static temperature [K]
    """

    T = Tt / (1.0 + 0.5 * (gamma - 1.0) * M**2)

    p = pt / (
        1.0 + 0.5 * (gamma - 1.0) * M**2
    ) ** (gamma / (gamma - 1.0))

    return p, T


def inlet_state(M):
    """
    Calculate the inlet static state and mass flow rate
    from the specified inlet stagnation conditions and Mach number.
    """

    p, T = static_from_total(M, pt_in, Tt_in)

    a = np.sqrt(gamma * R * T)
    u = M * a

    rho = p / (R * T)

    mdot = rho * u * A

    return p, T, a, u, rho, mdot


# =============================================================================
# RAYLEIGH RELATIONS
# =============================================================================

def rayleigh_pressure(M, C):
    """
    Rayleigh momentum relation:

        p (1 + gamma M^2) = constant = C

    """

    return C / (1.0 + gamma * M**2)


def rayleigh_temperature(M, C, mass_flux):
    """
    Calculate static temperature from the mass-conservation relation.

    G = rho u
      = p M sqrt(gamma / (R T))

    Therefore

        T = gamma p^2 M^2 / (R G^2)
    """

    p = rayleigh_pressure(M, C)

    T = (
        gamma
        * p**2
        * M**2
        / (R * mass_flux**2)
    )

    return T


def dT0_dM(M, T):
    """
    Derivative of total temperature with respect to Mach number.

    From the Rayleigh momentum and mass-conservation relations,

        dp/p = -2 gamma M / (1 + gamma M^2) dM

    and

        dT/T =
            2 [1/M - 2 gamma M/(1 + gamma M^2)] dM.

    Since

        T0 = T + u^2/(2 cp),

    the derivative is evaluated analytically.
    """

    # d(ln p)/dM
    dlnp_dM = (
        -2.0 * gamma * M
        / (1.0 + gamma * M**2)
    )

    # d(ln T)/dM
    dlnT_dM = 2.0 * (
        dlnp_dM + 1.0 / M
    )

    dT_dM = T * dlnT_dM

    # u^2 = gamma R T M^2
    u2 = gamma * R * T * M**2

    du2_dM = (
        gamma * R * T * M**2
        * (
            dlnT_dM
            + 2.0 / M
        )
    )

    # Tt = T + u^2/(2 cp)
    dTt_dM = (
        dT_dM
        + du2_dM / (2.0 * cp)
    )

    return dTt_dM


# =============================================================================
# INTEGRATE ONE RAYLEIGH CASE
# =============================================================================

def integrate_rayleigh(M_in, qw):
    """
    Integrate the Rayleigh flow from x=0 to x=L.

    State vector:
        y[0] = Mach number
        y[1] = total temperature

    The total temperature is included explicitly because

        cp dTt/dx = qw P / mdot.
    """

    # Inlet state
    p_in, T_in, a_in, u_in, rho_in, mdot = inlet_state(M_in)

    # Mass flux
    G = mdot / A

    # Rayleigh momentum constant
    C = p_in * (
        1.0 + gamma * M_in**2
    )

    # Heat addition per unit length
    qprime = qw * P

    def rhs(x, y):

        M = y[0]

        # Static pressure from momentum equation
        p = rayleigh_pressure(M, C)

        # Static temperature from mass conservation
        T = rayleigh_temperature(
            M,
            C,
            G
        )

        # dTt/dx from energy equation
        dTt_dx = qprime / (
            mdot * cp
        )

        # dTt/dM
        dTt_dM = dT0_dM(M, T)

        # Therefore
        #
        # dM/dx =
        #     (dTt/dx)/(dTt/dM)
        #
        dM_dx = dTt_dx / dTt_dM

        return np.array([
            dM_dx,
            dTt_dx
        ])

    solution = solve_ivp(
        rhs,
        [0.0, L],
        [M_in, Tt_in],
        rtol=RTOL,
        atol=ATOL,
        max_step=MAX_STEP
    )

    if not solution.success:
        raise RuntimeError(
            "Rayleigh integration failed: "
            + solution.message
        )

    return {
        "solution": solution,
        "mdot": mdot,
        "mass_flux": G,
        "momentum_constant": C,
        "p_in": p_in,
        "T_in": T_in,
        "a_in": a_in,
        "u_in": u_in,
        "rho_in": rho_in
    }


# =============================================================================
# EXTRACT COMPLETE FLOWFIELD
# =============================================================================

def postprocess(case):
    """
    Calculate all primitive and stagnation variables
    from the integrated Mach-number distribution.
    """

    sol = case["solution"]

    x = sol.t
    M = sol.y[0]
    Tt = sol.y[1]

    C = case["momentum_constant"]
    G = case["mass_flux"]

    # Static pressure
    p = rayleigh_pressure(M, C)

    # Static temperature
    T = rayleigh_temperature(
        M,
        C,
        G
    )

    # Speed of sound
    a = np.sqrt(gamma * R * T)

    # Velocity
    u = M * a

    # Density
    rho = p / (R * T)

    # Stagnation pressure is NOT constant in Rayleigh flow.
    #
    # It is calculated from the local static state and Mach number:
    #
    # pt = p (1 + (gamma-1)/2 M^2)^(gamma/(gamma-1))
    #
    pt = p * (
        1.0
        + 0.5 * (gamma - 1.0) * M**2
    ) ** (
        gamma / (gamma - 1.0)
    )

    # Total enthalpy
    ht = cp * Tt

    # Specific total enthalpy from static quantities
    ht_check = cp * T + 0.5 * u**2

    # Local heat flux
    qw = case["qw"]

    # Cumulative heat added per unit mass
    q_per_mass = (
        qw * P * x / case["mdot"]
    )

    return {
        "x": x,
        "M": M,
        "p": p,
        "T": T,
        "pt": pt,
        "Tt": Tt,
        "u": u,
        "a": a,
        "rho": rho,
        "ht": ht,
        "ht_check": ht_check,
        "q_per_mass": q_per_mass,
        "mdot": case["mdot"],
        "qw": qw
    }


# =============================================================================
# FIND INLET MACH NUMBER
# =============================================================================

def outlet_pressure_residual(M_in, qw):
    """
    Residual used to determine the inlet Mach number:

        F(M_in) = p(L) - p_out
    """

    # Zero heat flux has a particularly simple solution.
    if np.isclose(qw, 0.0):
        p_in, _, _, _, _, _ = inlet_state(M_in)
        return p_in - p_out

    case = integrate_rayleigh(
        M_in,
        qw
    )

    data = postprocess({
        **case,
        "qw": qw
    })

    return data["p"][-1] - p_out


def find_inlet_mach(qw):
    """
    Determine the subsonic inlet Mach number.

    A bracket search is performed first, followed by Brent's
    robust root-finding algorithm.
    """

    # For qw = 0 the solution is simply the isentropic
    # inlet condition satisfying p_in = p_out.
    if np.isclose(qw, 0.0):

        def residual(M):
            p, _, _, _, _, _ = inlet_state(M)
            return p - p_out

        return brentq(
            residual,
            M_MIN,
            M_MAX,
            xtol=1.0e-12,
            rtol=1.0e-12
        )

    # For heat addition, search for a subsonic root.
    M_grid = np.linspace(
        M_MIN,
        M_MAX,
        200
    )

    residuals = []

    for M in M_grid:
        try:
            residuals.append(
                outlet_pressure_residual(
                    M,
                    qw
                )
            )
        except Exception:
            residuals.append(np.nan)

    residuals = np.asarray(residuals)

    for i in range(len(M_grid) - 1):

        f1 = residuals[i]
        f2 = residuals[i + 1]

        if not np.isfinite(f1) or not np.isfinite(f2):
            continue

        if f1 * f2 < 0.0:

            M1 = M_grid[i]
            M2 = M_grid[i + 1]

            return brentq(
                lambda M: outlet_pressure_residual(
                    M,
                    qw
                ),
                M1,
                M2,
                xtol=1.0e-11,
                rtol=1.0e-11
            )

    raise RuntimeError(
        f"No subsonic solution found for "
        f"qw = {qw / 1e3:.3f} kW/m^2.\n"
        "The prescribed outlet pressure and heat flux "
        "may be incompatible with a subsonic Rayleigh-flow solution."
    )


# =============================================================================
# SOLVE ALL CASES
# =============================================================================

results = {}

for qw in qw_values:

    print()
    print("-" * 75)
    print(
        f"Solving qw = {qw / 1e3:.3f} kW/m^2"
    )
    print("-" * 75)

    # Determine inlet Mach number
    M_in = find_inlet_mach(qw)

    # Integrate
    case = integrate_rayleigh(
        M_in,
        qw
    )

    case["qw"] = qw
    case["M_in"] = M_in

    # Postprocess
    data = postprocess(case)

    results[qw] = data

    # Extract outlet quantities
    M_out = data["M"][-1]
    p_in = data["p"][0]
    p_final = data["p"][-1]

    T_in = data["T"][0]
    T_out = data["T"][-1]

    Tt_in = data["Tt"][0]
    Tt_out = data["Tt"][-1]

    pt_in_actual = data["pt"][0]
    pt_out = data["pt"][-1]

    mdot = data["mdot"]

    Q_total = qw * P * L

    print(f"Inlet Mach       : {M_in:.8f}")
    print(f"Outlet Mach      : {M_out:.8f}")
    print(f"Mass flow rate   : {mdot:.8f} kg/s")
    print()
    print(f"Inlet static p   : {p_in / 1e3:.6f} kPa")
    print(f"Outlet static p  : {p_final / 1e3:.6f} kPa")
    print()
    print(f"Inlet static T   : {T_in:.6f} K")
    print(f"Outlet static T  : {T_out:.6f} K")
    print()
    print(f"Inlet total T    : {Tt_in:.6f} K")
    print(f"Outlet total T   : {Tt_out:.6f} K")
    print()
    print(f"Inlet total p    : {pt_in_actual / 1e3:.6f} kPa")
    print(f"Outlet total p   : {pt_out / 1e3:.6f} kPa")
    print()
    print(f"Total heat added : {Q_total:.6f} W")


# =============================================================================
# SUMMARY TABLE
# =============================================================================

print()
print("=" * 100)
print("SUMMARY")
print("=" * 100)

print(
    f"{'qw [kW/m2]':>14}"
    f"{'M_in':>14}"
    f"{'M_out':>14}"
    f"{'mdot [kg/s]':>16}"
    f"{'Tt_out [K]':>16}"
    f"{'p_out [kPa]':>16}"
)

print("-" * 100)

for qw in qw_values:

    data = results[qw]

    print(
        f"{qw / 1e3:14.3f}"
        f"{data['M'][0]:14.8f}"
        f"{data['M'][-1]:14.8f}"
        f"{data['mdot']:16.8f}"
        f"{data['Tt'][-1]:16.6f}"
        f"{data['p'][-1] / 1e3:16.6f}"
    )


# =============================================================================
# PLOTS
# =============================================================================

# -------------------------------------------------------------------------
# Mach number
# -------------------------------------------------------------------------

plt.figure(figsize=(8, 5))

for qw in qw_values:

    data = results[qw]

    plt.plot(
        data["x"],
        data["M"],
        linewidth=2.0,
        label=(
            rf"$q_w={qw / 1e3:g}$ kW/m$^2$"
        )
    )

plt.xlabel("Axial coordinate, $x$ [m]")
plt.ylabel("Mach number, $M$")
plt.title("Rayleigh flow: Mach number")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()


# -------------------------------------------------------------------------
# Static pressure
# -------------------------------------------------------------------------

plt.figure(figsize=(8, 5))

for qw in qw_values:

    data = results[qw]

    plt.plot(
        data["x"],
        data["p"] / 1e3,
        linewidth=2.0,
        label=(
            rf"$q_w={qw / 1e3:g}$ kW/m$^2$"
        )
    )

plt.xlabel("Axial coordinate, $x$ [m]")
plt.ylabel("Static pressure, $p$ [kPa]")
plt.title("Rayleigh flow: static pressure")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()


# -------------------------------------------------------------------------
# Static temperature
# -------------------------------------------------------------------------

plt.figure(figsize=(8, 5))

for qw in qw_values:

    data = results[qw]

    plt.plot(
        data["x"],
        data["T"],
        linewidth=2.0,
        label=(
            rf"$q_w={qw / 1e3:g}$ kW/m$^2$"
        )
    )

plt.xlabel("Axial coordinate, $x$ [m]")
plt.ylabel("Static temperature, $T$ [K]")
plt.title("Rayleigh flow: static temperature")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()


# -------------------------------------------------------------------------
# Total temperature
# -------------------------------------------------------------------------

plt.figure(figsize=(8, 5))

for qw in qw_values:

    data = results[qw]

    plt.plot(
        data["x"],
        data["Tt"],
        linewidth=2.0,
        label=(
            rf"$q_w={qw / 1e3:g}$ kW/m$^2$"
        )
    )

plt.xlabel("Axial coordinate, $x$ [m]")
plt.ylabel("Total temperature, $T_t$ [K]")
plt.title("Rayleigh flow: total temperature")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()


# -------------------------------------------------------------------------
# Total pressure
# -------------------------------------------------------------------------

plt.figure(figsize=(8, 5))

for qw in qw_values:

    data = results[qw]

    plt.plot(
        data["x"],
        data["pt"] / 1e3,
        linewidth=2.0,
        label=(
            rf"$q_w={qw / 1e3:g}$ kW/m$^2$"
        )
    )

plt.xlabel("Axial coordinate, $x$ [m]")
plt.ylabel("Total pressure, $p_t$ [kPa]")
plt.title("Rayleigh flow: total pressure")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()


# -------------------------------------------------------------------------
# Velocity
# -------------------------------------------------------------------------

plt.figure(figsize=(8, 5))

for qw in qw_values:

    data = results[qw]

    plt.plot(
        data["x"],
        data["u"],
        linewidth=2.0,
        label=(
            rf"$q_w={qw / 1e3:g}$ kW/m$^2$"
        )
    )

plt.xlabel("Axial coordinate, $x$ [m]")
plt.ylabel("Velocity, $u$ [m/s]")
plt.title("Rayleigh flow: velocity")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

# save all results in a single pickle, containting a dictionary with the results for each qw
import pickle
with open('analytical_results.pkl', 'wb') as f:
    pickle.dump(results, f)


plt.show()